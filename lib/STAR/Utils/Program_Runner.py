import os.path
import subprocess


class Program_Runner:

    def __init__(self, cmd, scratch_dir):
        self.scratch_dir = scratch_dir
        self.executableName = cmd

    def run(self, command, cwd_dir=None):
        ''' options is an array of command-line parameters passed to the RQCFilter App '''
        cmmd = command

        if not cwd_dir:
          cwd_dir = self.scratch_dir

        print('\nRunning: ' + ' '.join(cmmd))
        p = subprocess.Popen(cmmd, cwd=cwd_dir, shell=False)
        exitCode = p.wait()

        if (exitCode == 0):
            print('\n' + ' '.join(cmmd) + ' was executed successfully, exit code was: ' + str(exitCode))
        else:
            raise ValueError('Error running command: ' + ' '.join(cmmd) + '\n' +
                             'Exit Code: ' + str(exitCode))

        return exitCode

