#!/usr/bin/env python3

import sys

tasks = {
    'scheme_create': 'Create a bioHansel compatible genotyping scheme',
    'scheme_evaluate': 'Evaluate genotyping scheme results for a test panel',
    'scheme_validate': 'Test genotyping scheme using bioHansel',
    'test': 'Test bioCanon functionality on a small dataset',
    'version': 'Print version and exit',
}


ordered_tasks = [
    'scheme_create',
    'scheme_evaluate',
    'scheme_validate',
    'test',
    'version'
]


def print_usage_and_exit():
    print('Usage: bioCanon <command> [options] <required arguments>', file=sys.stderr)
    print('\nTo get minimal usage for a command use:\nbioCanon command', file=sys.stderr)
    print('\nTo get full help for a command use one of:\nbioCanon command -h\nbioCanon command --help\n', file=sys.stderr)
    print('\nAvailable commands:\n', file=sys.stderr)
    max_task_length = max([len(x) for x in list(tasks.keys())]) + 1
    for task in ordered_tasks:
        print('{{0: <{}}}'.format(max_task_length).format(task), tasks[task], sep=' ', file=sys.stderr)
    sys.exit(0)

if len(sys.argv) == 1 or sys.argv[1] in ['-h', '-help', '--help']:
    print_usage_and_exit()

task = sys.argv.pop(1)

if task not in tasks:
    print('Task "' + task + '" not recognised. Cannot continue.\n', file=sys.stderr)
    print_usage_and_exit()



exec('import bioCanon.' + task)
exec('bioCanon.' + task + '.run()')