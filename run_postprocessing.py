"""使用共享后处理分析已保存的 WEIDA 结果。"""
import argparse
from pathlib import Path

from weida_postprocessing import save_results
from weida_postprocessing.adapters import load_weida


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--features', type=Path, help='可选的 run_features 输出目录')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    table = load_weida(args.input, args.features)
    print(save_results(table, args.output, input_path=args.input, features_path=args.features))


if __name__ == '__main__':
    main()
