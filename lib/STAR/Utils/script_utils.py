import logging
import os
import subprocess
import traceback
from zipfile import ZipFile
from os import listdir
from os.path import isfile, join


'''
A utility python module containing a set of methods necessary for this kbase
module.
'''

LEVELS = {'debug': logging.DEBUG,
          'info': logging.INFO,
          'warning': logging.WARNING,
          'error': logging.ERROR,
          'critical': logging.CRITICAL}


def create_logger(log_dir, name):
    """Create a logger
    args: name (str): name of logger
    returns: logger (obj): logging.Logger instance
    """
    logger = logging.getLogger(name)
    fmt = logging.Formatter('%(asctime)s - %(process)d - %(name)s - '
                            ' %(levelname)s -%(message)s')
    hdl = logging.FileHandler(os.path.join(log_dir, name + '.log'))
    hdl.setFormatter(fmt)

    logger.addHandler(hdl)

    return logger


def if_obj_exists(logger, ws_client, ws_id, o_type, obj_l):
    obj_list = ws_client.list_objects({"workspaces": [ws_id], "type": o_type, 'showHidden': 1})
    obj_names = [i[1] for i in obj_list]
    existing_names = [i for i in obj_l if i in obj_names]
    obj_ids = None
    if len(existing_names) != 0:
        e_queries = [{'name': j, 'workspace': ws_id} for j in existing_names]
        e_infos = ws_client.get_object_info_new({"objects": e_queries})
        obj_ids = [(str(k[1]), (str(k[6]) + '/' + str(k[0]) + '/' + str(k[4]))) for k in e_infos]
    return obj_ids


def log(message, level=logging.INFO, logger=None):
    if logger is None:
        if level == logging.DEBUG:
            print('\nDEBUG: ' + message + '\n')
        elif level == logging.INFO:
            print('\nINFO: ' + message + '\n')
        elif level == logging.WARNING:
            print('\nWARNING: ' + message + '\n')
        elif level == logging.ERROR:
            print('\nERROR: ' + message + '\n')
        elif level == logging.CRITICAL:
            print('\nCRITICAL: ' + message + '\n')
    else:
        logger.log(level, '\n' + message + '\n')


def zip_files(logger, src_path, output_fn):
    """
    Compress all index files (not directory) into an output zip file on disk.
    """

    files = [f for f in listdir(src_path) if isfile(join(src_path, f))]
    with ZipFile(output_fn, 'w', allowZip64=True) as izip:
        for f in files:
            izip.write(join(src_path, f), f)


def unzip_files(logger, src_fn, dst_path):
    """
    Extract all index files into an output zip file on disk.
    """

    with ZipFile(src_fn, 'r') as ozip:
        ozip.extractall(dst_path)


def whereis(program):
    """
    returns path of program if it exists in your ``$PATH`` variable or `
    `None`` otherwise
    """
    for path in os.environ.get('PATH', '').split(':'):
        if os.path.exists(os.path.join(path, program)) and not os.path.isdir(
                os.path.join(path, program)):
            return os.path.join(path, program)
    return None


def runProgram(logger=None,
               progName=None,
               argStr=None,
               script_dir=None,
               working_dir=None):
    """
    Convenience func to handle calling and monitoring output of external programs.
    :param progName: name of system program command
    :param argStr: string containing command line options for ``progName``
    :returns: subprocess.communicate object
    """

    # Ensure program is callable.
    if script_dir is not None:
        progPath = os.path.join(script_dir, progName)
    else:
        progPath = progName
    progPath = whereis(progName)
    if not progPath:
        raise RuntimeError(
            None,
            '{0} command not found in your PATH environmental variable. {1}'.format(
                progName,
                os.environ.get(
                    'PATH',
                    '')))

    # Construct shell command
    cmdStr = "%s %s" % (progPath, argStr)
    print "Executing : " + cmdStr
    if logger is not None:
        logger.info("Executing : " + cmdStr)
    # if working_dir is None:
        logger.info("Executing: " + cmdStr + " on cwd")
    else:
        logger.info("Executing: " + cmdStr + " on " + working_dir)

    # Set up process obj
    process = subprocess.Popen(cmdStr,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=working_dir)
    # Get results
    result, stderr = process.communicate()
    # print result
    # print stderr
    # keep this until your code is stable for easier debugging
    if logger is not None and result is not None and len(result) > 0:
        logger.info(result)
    else:
        print result
    if logger is not None and stderr is not None and len(stderr) > 0:
        logger.info(stderr)
    else:
        print stderr

    # Check returncode for success/failure
    if process.returncode != 0:
        raise Exception("Command execution failed  {0}".format(
            "".join(traceback.format_exc())))
        raise RuntimeError(
            'Return Code : {0} , result {1} , progName {2}'.format(
                process.returncode, result, progName))

    # Return result
    return {"result": result, "stderr": stderr}


def check_sys_stat(logger):
    check_disk_space(logger)
    check_memory_usage(logger)
    check_cpu_usage(logger)


def check_disk_space(logger):
    runProgram(logger=logger, progName="df", argStr="-h")


def check_memory_usage(logger):
    runProgram(logger=logger, progName="vmstat", argStr="-s")


def check_cpu_usage(logger):
    runProgram(logger=logger, progName="mpstat", argStr="-P ALL")

