'''
Usage: python extract_function.py --from A.ipynb --to A.py B.py C.py
'''
if __name__ == '__main__':
    import argparse, re
    import json
    parser = argparse.ArgumentParser()
    parser.add_argument("--from", dest = "source")
    parser.add_argument("--to", nargs = '+')
    args = parser.parse_args()
    print("Extracting from {}".format(args.source))
    output = []
    source = []
    harvest = False
    for cell in json.load(open(args.source))['cells']:
        if cell['cell_type'] == 'code':
            source.extend([x.rstrip() for x in cell['source']])
    for line in source:
        if line.startswith('#'):
            continue
        if line.startswith('import') or line.startswith('from'):
            output.append(line)
        if line.startswith('def') or line.startswith('class'):
            harvest = True
            output.append(line)
            continue
        if harvest:
            if line.lstrip() != line:
                output.append(line)
            else:
                output.append('\n')
                harvest = False
    for item in args.to:
        with open(item, 'w') as f:
            f.write('\n'.join(output))
    print("Done!")
