import os

def create_dirs(args):
    os.mkdir(f"work/{args.project}", exists=True)
    os.mkdir(f"results/{args.project}", exists=True)