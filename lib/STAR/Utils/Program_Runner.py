import os.path
import subprocess


class Program_Runner:


    def __init__(self, cmd, scratch_dir):
        self.scratch_dir = scratch_dir
        self.executableName = cmd

    def run(self, command, cwd=None):
        ''' options is an array of command-line parameters passed to the RQCFilter App '''
        command = [os.path.join(self.executableName, command)]

        print('In working directory: ' + ' '.join(command))
        print('Running: ' + ' '.join(command))


        if not cwd:
          cwd = self.scratch_dir

        p = subprocess.Popen(command, cwd=cwd, shell=False)
        exitCode = p.wait()

        if (exitCode == 0):
            print(command + ' was executed successfully, exit code was: ' + str(exitCode))
        else:
            raise ValueError('Error running command: ' + ' '.join(command) + '\n' +
                             'Exit Code: ' + str(exitCode))

        return exitCode

