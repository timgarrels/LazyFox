import argparse
import os

from notebooks.Datasets import download, urls
from os import makedirs


def main():
    parser = argparse.ArgumentParser(description='Download and Extract SNAP Datasets')
    parser.add_argument('--dataset', help='The dataset to download', required=True, choices=urls.keys())
    parser.add_argument('--output', help='The directory to download to', required=True)

    args = parser.parse_args()

    if args.dataset not in urls.keys():
        exit(f'No such dataset {args.dataset}')

    os.makedirs(args.output, exist_ok=True)

    download(args.dataset, args.output)


if __name__ == '__main__':
    main()
