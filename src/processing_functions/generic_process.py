import os

def create_dirs(args):
    os.makedirs(f"work/{args.project}", exist_ok=True)
    os.makedirs(f"work/{args.project}/samplesheets", exist_ok=True)
    os.makedirs(f"results/{args.project}", exist_ok=True)
