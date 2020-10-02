import os
import time
import matplotlib.pyplot as plt


EH_RLOOPER = '/home/ethan/Documents/github/rlooper/bin/rlooper'
OG_RLOOPER = '/home/ethan/Documents/github/OG_rlooper/rlooper/bin/rlooper'
OUTPUT_FOLDER = '/home/ethan/Documents/github/rlooper/testing/time_trials/time_trial_output'
FASTA = '/home/ethan/Documents/github/rlooper/testing/time_trials/test.fasta'

if not os.path.exists(OUTPUT_FOLDER): 
    os.mkdir(OUTPUT_FOLDER)

SEQ_HEADER = ">HG19_AIRN_PFC53_REVERSE_dna range=PFC53FIXED:1-{} 5'pad=0 3'pad=0 strand=- repeatMasking=none"
length_log, og_runtime, eh_runtime = [], [], []


def run_looper(exe_path, iden):
    
    output = os.path.join(OUTPUT_FOLDER, iden)
    if not os.path.isdir(output):
        os.mkdir(output)
    
    cmd = '{} {} {} --sigma -0.10 --dump'.format(exe_path, FASTA, output)
    print('CMD:', cmd)
    s = time.time()
    status = os.system(cmd)
    e = time.time()
    print(iden, status)

    return e - s # time to run in seconds


def make_test_seq(seq_len):
    header = SEQ_HEADER.format(seq_len)
    os.system('/home/ethan/bin/random_dna {} --header "{}" > test.fasta'.format(seq_len, header))


def time_trial(n=10):
    for i in range(n):
        seq_len = 500 * (i+1)  # found bug where program stuck in loop is seq len = 1
        print('RUNNING FOR SEQ LENGTH', seq_len, '\n\n')
        make_test_seq(seq_len)
        print('Made fresh test sequence')
        length_log.append(seq_len)
        eh_runtime.append(run_looper(EH_RLOOPER, "eh_test"))
        og_runtime.append(run_looper(OG_RLOOPER, "og_test"))

def plot_results(lens, og_times, eh_times):
    plt.plot(lens, eh_times, label='Dynamic Free Energy Calculations')
    plt.plot(lens, og_times, label='Stock RLooper')
    plt.legend(loc="upper right")

    plt.xlabel('Sequence length (bp)')
    plt.ylabel('Runtime in seconds')
    #plt.show()

time_trial(8)
plot_results(length_log, og_runtime, eh_runtime)
plt.savefig('/home/ethan/Documents/github/lab_notes/Chedin/rlooper_optimization/dynamic_energy_calcs/compare_2.png')




