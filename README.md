您可以将基因序列的对齐文件放置于文件夹中，指定文件夹运行，但同时您必须保证您的文件的扩展名为fa/fas/fasta，您也可以只输入一个序列比对文件，在gene模式下
这个脚本会计算每一个位点的碱基频率，如果您指定了-d，他会计算核苷酸的多样性（Pi）值，记得使用-t来指定输出文件。当然，对于每一个文件您也可以使用滑窗方法
指定步长（-s）来拆分序列，以计算每一个划窗中的碱基频率和Pi值，这里只需要使用-m window来指定模式和-s 来指定步长。


python /path/script/CalculationPi_seq.py -f ./alignSeq/ -o ./SeqPi/ -c acgt -d -t Seqpi.tsv -m gene


python /path/script/CalculationPi_seq.py -f ./concatcp.fasta -o ./SeqPi/ -c acgt -d -t Seqpi02.tsv -m window -s 1000


请注意，该版本是在python3.12来编写和测试的，运行前请检查您的python版本，一般来说python3都可以运行，因为这个脚本比较简单，如果报错请检查您的文件和相关
输入指令是否正确，没有基金支持该项目进行维护，所以基本情况下都不会维护和回复问题，当然如果您发现bug恳请您留言！

如果您能引用就太好啦！

祝您工作和生活顺利，科研生活精彩！！！


FASTA 位点频率与多样性计算工具，支持 gene 与 window 模式


options:

  -h, --help            show this help message and exit
  
  -f FASTA, --fasta FASTA
  
                        输入FASTA文件或目录路径
                        
  -o OUTDIR, --outdir OUTDIR
  
                        输出目录 [默认: ./out]
                        
  -c CODE, --code CODE  计数的碱基字符，默认: ACGT
  
  -g GAPCODE, --gapcode GAPCODE
  
                        gap字符，默认: '-'
                        
  -d, --div             是否计算多样性 Pi 值
  
  -t DIVOUT, --divout DIVOUT
  
                        多样性输出文件路径
                        
  -m {gene,window}, --model {gene,window}
  
                        运行模式：gene（目录模式）或 window（单文件滑窗）
                        
  -s SIZE, --size SIZE  滑窗大小，仅 window 模式生效，默认: 500
