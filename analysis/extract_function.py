import argparse, re
'''
Usage: python extract_function.py --from A.py --to B.py
'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--from", dest = "source")
    parser.add_argument("--to")
    args = parser.parse_args()
    output = []
    harvest = False
    for line in open(args.source).readlines():
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
    with open(args.to, 'w') as f:
        f.write(''.join(output))
