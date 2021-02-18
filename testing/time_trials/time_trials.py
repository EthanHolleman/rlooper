import os
import time


EH_RLOOPER = ''
OG_RLOOPER = ''
OUTPUT_FOLDER = 'time_trial_output'
FASTA = os.path.join(OUTPUT_FOLDER, 'test.fasta')

os.mkdir(OUTPUT_FOLDER)
SEQ_HEADER = ">HG19_AIRN_PFC53_REVERSE_dna range=PFC53FIXED:1-{} 5'pad=0 3'pad=0 strand=- repeatMasking=none"
length_log, og_runtime, eh_runtime = [], [], []


def run_looper(exe_path, iden):
    s = time.time()
    output = os.path.join(OUTPUT_FOLDER, iden)
    if not os.path.isdir(output):
        os.mkdir(output)
    
    cmd = os.system(
        '{} {}'
    )
    e = time.time()




for i in range(0, 10):
    n = 500 * i
    # prepare the DNA sequence for testing
    header = SEQ_HEADER.format(n)
    os.system('random_dna {} --header {} > test.fasta'.format(n, header))

    # run OG program

    start = time.time()
    os.system('{} {} {} --sigma -0.10'.format(OG_RLOOPER, FASTA, ))

    test.fasta rlooper_test_dynamic_output --sigma -0.10


