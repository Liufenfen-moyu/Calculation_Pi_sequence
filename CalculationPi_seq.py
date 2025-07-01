import os
import argparse
import numpy as np
from numba import njit

def make_fastagrid(fasta_path):
    """
    将FASTA文件读取为 numpy array，每行为一个序列，每列为一个碱基位置。
    """
    seqs = []
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_seq:
                    seqs.append("".join(current_seq).upper())
                    current_seq = []
            else:
                current_seq.append(line)
        if current_seq:
            seqs.append("".join(current_seq).upper())

    if len(set(len(s) for s in seqs)) != 1:
        raise ValueError(f"FASTA 文件中序列长度不一致: {fasta_path}")

    # 转换为numpy二维矩阵
    seq_array = np.array([list(seq) for seq in seqs])
    return seq_array


def score_fastagrid(grid, code):
    """
    对每个位点统计各碱基的出现次数，以及除目标碱基以外的字符数。
    返回结果为 numpy array，每行：[位置, A数, C数, G数, T数, ..., Other数]。
    """
    score_grid = []
    for idx, column in enumerate(grid.T, start=1):  # 转置得到列为位点
        counts = [idx] + [np.sum(column == base) for base in code]
        counts.append(len(column) - sum(counts[1:]))  # 统计非指定字符数量
        score_grid.append(counts)
    return np.array(score_grid)


@njit
def calculate_X(score_grid):
    """
    使用 Numba 加速的多样性统计量计算函数。
    score_grid: numpy 2D array，每行：[位置, A数, C数, G数, T数, ..., Other数]
    返回值：
    - X: 所有位点的 minor allele 次数之和
    - L: 比对长度（即位点数）
    - D: 序列数（= 每个位点的计数和）
    """
    L = score_grid.shape[0]
    D = 0
    X = 0.0

    for i in range(L):
        row = score_grid[i]
        allele_counts = row[1:-1]  # 去除位置与Other列
        row_sum = 0
        max_val = 0
        for val in allele_counts:
            row_sum += val
            if val > max_val:
                max_val = val
        D = row_sum  # 假设每个位点都有相同的序列数
        minor = row_sum - max_val
        X += minor

    return X, float(L), float(D)


def run_all(options):
    fastadir = options.fastadir
    outdir = options.outdir
    code = list(options.code.upper())
    gapcode = options.gapcode
    code_with_gap = code + [gapcode]
    calculate_div = options.div
    divout_file = options.divout

    print(f"输入FASTA目录: {fastadir}")
    print(f"输出目录: {outdir}")
    print(f"计数字符: {code}（含gap符号: {gapcode}）")
    print(f"是否计算多样性: {calculate_div}")

    if not os.path.isdir(fastadir):
        raise FileNotFoundError(f"找不到FASTA目录: {fastadir}")

    fasta_files = [f for f in os.listdir(fastadir) if f.lower().endswith(('.fa', '.fasta'))]
    if not fasta_files:
        raise FileNotFoundError("FASTA目录中未找到.fasta或.fa文件")

    os.makedirs(outdir, exist_ok=True)

    diversity_output = []
    if calculate_div:
        diversity_output.append("File\tX\tLength\tDepth\tPi")

    for fasta in fasta_files:
        fasta_path = os.path.join(fastadir, fasta)
        try:
            grid = make_fastagrid(fasta_path)
        except Exception as e:
            print(f"[跳过] 读取失败 {fasta_path}：{e}")
            continue

        score_grid = score_fastagrid(grid, code_with_gap)

        # 保存位点频率表
        out_path = os.path.join(outdir, f"{fasta}.txt")
        header = "Position\t" + "\t".join(code_with_gap) + "\tOther"
        np.savetxt(out_path, score_grid, fmt="%i", delimiter="\t", header=header, comments="")

        # 多样性计算
        if calculate_div:
            try:
                X, L, D = calculate_X(score_grid)
                pi = X / (((D * (D - 1)) / 2) * L) if D > 1 and L > 0 else 0.0
                diversity_output.append(f"{fasta}\t{int(X)}\t{int(L)}\t{int(D)}\t{pi:.6f}")
            except Exception as e:
                print(f"[警告] 多样性计算失败: {fasta}: {e}")

    # 输出多样性文件
    if calculate_div:
        with open(divout_file, "w") as f:
            f.write("\n".join(diversity_output))
        print(f"多样性输出写入: {divout_file}")


def main():
    parser = argparse.ArgumentParser(description="统计FASTA文件中每个位点的碱基频率，并可选计算碱基多样性指标。")
    parser.add_argument("-f", "--fastadir", required=True,
                        help="包含FASTA文件的输入目录")
    parser.add_argument("-o", "--outdir", default="out",
                        help="输出目录，默认为 ./out")
    parser.add_argument("-c", "--code", default="ACGT",
                        help="统计的碱基代码，默认为 'ACGT'")
    parser.add_argument("-g", "--gapcode", default="-",
                        help="表示gap的字符，默认 '-'")
    parser.add_argument("-d", "--div", action="store_true",
                        help="是否计算碱基多样性 Pi 值")
    parser.add_argument("-t", "--divout", default="diversity.tsv",
                        help="输出的多样性文件路径，仅在指定 -d 时有效")

    args = parser.parse_args()
    run_all(args)


if __name__ == "__main__":
    main()

