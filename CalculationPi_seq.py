import os
import argparse
import numpy as np
from numba import njit

@njit
def calculate_X(score_grid):
    L = score_grid.shape[0]
    D = 0
    X = 0.0
    for i in range(L):
        row = score_grid[i]
        allele_counts = row[1:-1]
        row_sum = 0
        max_val = 0
        for val in allele_counts:
            row_sum += val
            if val > max_val:
                max_val = val
        D = row_sum
        minor = row_sum - max_val
        X += minor
    return X, float(L), float(D)

def make_fastagrid(fasta_path):
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
    return np.array([list(seq) for seq in seqs])

def score_fastagrid(grid, code):
    score_grid = []
    for idx, column in enumerate(grid.T, start=1):
        counts = [idx] + [np.sum(column == base) for base in code]
        counts.append(len(column) - sum(counts[1:]))
        score_grid.append(counts)
    return np.array(score_grid)

def write_score_and_diversity(score_grid, code, output_path, fasta_name, calc_div=False, diversity_output=None):
    header = "Position\t" + "\t".join(code) + "\tOther"
    np.savetxt(output_path, score_grid, fmt="%i", delimiter="\t", header=header, comments="")
    if calc_div:
        X, L, D = calculate_X(score_grid)
        pi = X / (((D * (D - 1)) / 2) * L) if D > 1 and L > 0 else 0.0
        diversity_output.append(f"{fasta_name}\t{int(X)}\t{int(L)}\t{int(D)}\t{pi:.6f}")

def run_gene_model(args, code_with_gap):
    fasta_dir = args.fasta
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    diversity_output = ["File\tX\tLength\tDepth\tPi"] if args.div else []
    fasta_files = [f for f in os.listdir(fasta_dir) if f.lower().endswith(('.fa', '.fasta', 'fas'))]
    if not fasta_files:
        raise FileNotFoundError("未找到FASTA文件")
    for fasta in fasta_files:
        fasta_path = os.path.join(fasta_dir, fasta)
        try:
            grid = make_fastagrid(fasta_path)
            score_grid = score_fastagrid(grid, code_with_gap)
            out_file = os.path.join(outdir, f"{fasta}.txt")
            write_score_and_diversity(score_grid, code_with_gap, out_file, fasta, args.div, diversity_output)
        except Exception as e:
            print(f"[跳过] 文件 {fasta}: {e}")
    if args.div:
        with open(args.divout, "w") as f:
            f.write("\n".join(diversity_output))
        print(f"多样性输出：{args.divout}")

def run_window_model(args, code_with_gap):
    fasta_file = args.fasta
    outdir = args.outdir
    win_size = args.size
    os.makedirs(outdir, exist_ok=True)
    grid = make_fastagrid(fasta_file)
    n_seqs, seq_len = grid.shape
    diversity_output = ["File\tX\tLength\tDepth\tPi"] if args.div else []
    fasta_name = os.path.basename(fasta_file)

    # 滑窗切片
    for start in range(0, seq_len, win_size):
        end = min(start + win_size, seq_len)
        if end - start < 2:
            continue
        window_grid = grid[:, start:end]
        score_grid = score_fastagrid(window_grid, code_with_gap)
        window_name = f"{fasta_name}_win_{start+1}_{end}"
        out_file = os.path.join(outdir, f"{window_name}.txt")
        write_score_and_diversity(score_grid, code_with_gap, out_file, window_name, args.div, diversity_output)

    if args.div:
        with open(args.divout, "w") as f:
            f.write("\n".join(diversity_output))
        print(f"多样性输出：{args.divout}")

def main():
    parser = argparse.ArgumentParser(description="FASTA 位点频率与多样性计算工具，支持 gene 与 window 模式")
    parser.add_argument("-f", "--fasta", required=True,
                        help="输入FASTA文件或目录路径")
    parser.add_argument("-o", "--outdir", default="out",
                        help="输出目录 [默认: ./out]")
    parser.add_argument("-c", "--code", default="ACGT",
                        help="计数的碱基字符，默认: ACGT")
    parser.add_argument("-g", "--gapcode", default="-",
                        help="gap字符，默认: '-'")
    parser.add_argument("-d", "--div", action="store_true",
                        help="是否计算多样性 Pi 值")
    parser.add_argument("-t", "--divout", default="diversity.tsv",
                        help="多样性输出文件路径")
    parser.add_argument("-m", "--model", choices=["gene", "window"], default="gene",
                        help="运行模式：gene（目录模式）或 window（单文件滑窗）")
    parser.add_argument("-s", "--size", type=int, default=500,
                        help="滑窗大小，仅 window 模式生效，默认: 500")

    args = parser.parse_args()
    code_with_gap = list(args.code.upper()) + [args.gapcode]

    if args.model == "gene":
        if not os.path.isdir(args.fasta):
            raise ValueError("gene 模式下输入必须为目录")
        run_gene_model(args, code_with_gap)
    elif args.model == "window":
        if not os.path.isfile(args.fasta):
            raise ValueError("window 模式下输入必须为单个FASTA文件")
        run_window_model(args, code_with_gap)

if __name__ == "__main__":
    main()
