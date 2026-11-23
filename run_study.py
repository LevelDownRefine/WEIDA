"""统一运行 WEiDA / PhiID：读取 → 幅度预处理 → 算法 → 三组比较 → 保存与作图。"""
import argparse
from research_core.study import run_study


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('--method', required=True, choices=['weida', 'phiid'])
    parser.add_argument('--preprocessing', required=True, choices=['loaded', 'weida_amplitude'])
    parser.add_argument('--tr', type=float, default=2.1)
    parser.add_argument('--exclude-stn', nargs='*', default=[])
    parser.add_argument('--expected-regions', type=int, default=78)
    parser.add_argument('--no-plots', action='store_true')
    parser.add_argument('--k', type=int, help='WEiDA only; default 4')
    parser.add_argument('--seed', type=int, help='WEiDA only; default 0')
    parser.add_argument('--replicates', type=int, help='WEiDA only; default 300')
    parser.add_argument('--lag', type=int, help='PhiID only; default 1 frame')
    args = parser.parse_args()
    options = {name: getattr(args, name) for name in ('k', 'seed', 'replicates', 'lag')
               if getattr(args, name) is not None}
    allowed = {'weida': {'k', 'seed', 'replicates'}, 'phiid': {'lag'}}
    if set(options) - allowed[args.method]:
        parser.error(f'参数不适用于 {args.method}: {sorted(set(options) - allowed[args.method])}')
    run_study(args.input, args.output, method=args.method, preprocessing=args.preprocessing,
              tr=args.tr, exclude_stn=args.exclude_stn, expected_regions=args.expected_regions,
              options=options, plots=not args.no_plots)
    print(f'{args.method}: {args.output}')


if __name__ == '__main__':
    main()
