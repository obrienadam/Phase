import os
import re
import argparse
import tarfile
import time

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('ARCHIVE', type=str)
    parser.add_argument('--dir', type=str, help='directory containing the solution files')
    parser.add_argument('--delay', type=float, default=3., help='delay in checking for new files')

    args = parser.parse_args()

    prog = re.compile(r'^Proc\d+$')
    solution_set = set()

    with tarfile.open(args.ARCHIVE, 'w') as f:
        for dir in (dir for dir in os.listdir(args.dir) if prog.match(dir) and os.path.isdir(os.path.join(args.dir, dir))):
            f.add(os.path.join(args.dir, dir))
            solution_set.add(dir)

    prog = re.compile(r'^\d+(?:\.\d+)?$')

    while True:
        new_solns = set(f for f in os.listdir(args.dir) if
                        prog.match(f) and os.path.isdir(os.path.join(args.dir, f))) - solution_set

        if len(new_solns) >= 2:
            new_solns = list(new_solns)
            new_solns.sort(key=lambda x: float(x))

            # Never use the last file in case it is still being written
            new_solns.pop(-1)

            with tarfile.open(args.ARCHIVE, 'a') as f:
                for dir in new_solns:
                    f.add(os.path.join(args.dir, dir))

            solution_set.update(new_solns)

            print('{} new files added to {}.'.format(len(new_solns), args.ARCHIVE))
        else:
            print('No new files added to {}.'.format(args.ARCHIVE))

        print('Sleeping for {} seconds...'.format(args.delay))
        time.sleep(args.delay)
