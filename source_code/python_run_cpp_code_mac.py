import subprocess
#import sys

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)
def run_spectrometer_code():
    run_args = ("./run")
    #Or just:
    #args = "bin/bar -c somefile.xml -d text.txt -r aString -f anotherString".split()
    popen = subprocess.Popen(run_args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    print(output)