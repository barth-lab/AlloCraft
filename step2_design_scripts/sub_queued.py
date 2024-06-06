
import re
import os
import subprocess

from pathlib import Path


queue_file_path = Path('/path/to/ism/cluster_queue.txt')
queue_file_data = queue_file_path.read_text().split('\n')

new_processes   = 0
blank_lines     = 0


for line in queue_file_data:

    if line == '':
        blank_lines += 1
        continue

    if new_processes >= 1000:
    	break
    
    command  = 'sbatch {}'.format(line)

    print(command)

    command_process = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
    print(command_process.stdout.decode().strip(), )

    new_processes += 1


print('New processes added:', new_processes)


queue_file = open(queue_file_path, mode='w+')

for line in queue_file_data[(new_processes+blank_lines):]:
    queue_file.write(line)
    queue_file.write('\n')

queue_file.close()

