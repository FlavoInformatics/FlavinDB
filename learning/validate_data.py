import sys
import pandas as pd
from redox_potential_helpers import trim_list

def main():
    raise Warning("Checking data is still in the initial stages and is definitely not totally bulletproof - please add to it as needed.")
    if len(sys.argv) < 3:
        print("./test.py filename category")
    filename = sys.argv[1]
    category = sys.argv[2]

    names = list(pd.read_csv(filename)[category].unique())

    found, not_found = trim_list(names, category)
    print('(unique) total count:', len(names))
    print('\tfound:', len(found))
    print('\tnot_found:', len(not_found))

    for name in not_found:
        print('\t', name)


if __name__ == '__main__':
    main()
