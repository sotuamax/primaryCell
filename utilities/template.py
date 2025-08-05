#!/usr/bin/env python3

import argparse 

def args_parser():
    parser = argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, add_help = False)
    parser.add_argument("input", help = "input file")
    parser.add_argument("-o", "--output", help = "output file prefix name")
    args = parser.parse_args()
    return args

def main():
    args = args_parser()
    args.input 
    args.output 

if __name__ == "__main__":
    main()

