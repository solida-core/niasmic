import errno
import multiprocessing
import os.path
import psutil


def _gatk_multi_arg(flag, args):
    flag += " "
    return " ".join(flag + arg for arg in args)

def _multi_flag(arguments):
    if arguments:
        return " ".join(flag + " " + arg for flag, arg in arguments)
    return ''

def _get_samples_set(samples_files, flag='-sf'):
    arguments = []
    set_arg=config.get("samples_set", None)
    if set_arg and set_arg in samples_files:
        return "".join(flag + " " + samples_files[set_arg])
    for samples_set in samples_files.keys():
        arguments.append("".join(flag + " " + samples_files[samples_set]))
    return ' '.join(arguments)

def total_physical_mem_size():
    mem = psutil.virtual_memory()
    return mem.total

def cpu_count():
    return multiprocessing.cpu_count()


def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT)+" (path must be absolute)", filepath)
    return filepath

def conservative_cpu_count(reserve_cores=1, max_cores=5):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)

def tmp_path(path=''):
    """
    if does not exists, create path and return it. If any errors, return
    default path
    :param path: path
    :return: path
    """
    default_path = os.getenv('TMPDIR', '/tmp')
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                return default_path
        return path
    return default_path

def java_params(tmp_dir='', percentage_to_preserve=20, stock_mem=1024 ** 3,
                stock_cpu=2, multiply_by=1):
    """
    Set Java params
    :param tmp_dir: path to tmpdir
    :param percentage_to_preserve: percentage of resources to preserve
    :param stock_mem: min memory to preserve
    :param stock_cpu: min cpu to preserve
    :param multiply_by: multiply base resource by this param
    :return: string to return to configure java environments
    """

    def bytes2human(n):
        # http://code.activestate.com/recipes/578019
        # >>> bytes2human(10000)
        # '9.8K'
        # >>> bytes2human(100001221)
        # '95.4M'
        symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
        prefix = {}
        for i, s in enumerate(symbols):
            prefix[s] = 1 << (i + 1) * 10
        for s in reversed(symbols):
            if n >= prefix[s]:
                value = float(n) / prefix[s]
                return '%.0f%s' % (value, s)
        return "%sB" % n

    def preserve(resource, percentage, stock):
        preserved = resource - max(resource * percentage // 100, stock)
        return preserved if preserved != 0 else stock

    # def preserve(resource, percentage, stock):
    #     return resource - max(resource * percentage // 100, stock)

    params_template = "'-Xms{} -Xmx{} -XX:ParallelGCThreads={} " \
                      "-Djava.io.tmpdir={}'"

    mem_min = 1024 ** 3 * 2  # 2GB

    mem_size = preserve(total_physical_mem_size(), percentage_to_preserve,
                        stock_mem)
    cpu_nums = preserve(cpu_count(), percentage_to_preserve, stock_cpu)
    tmpdir = tmp_path(tmp_dir)

    return params_template.format(bytes2human(mem_min).lower(),
                                  bytes2human(max(mem_size//cpu_nums*multiply_by,
                                                  mem_min)).lower(),
                                  min(cpu_nums, multiply_by),
                                  tmpdir)


# def references_abs_path(ref='references'):
#     references = config.get(ref)
#     basepath = references['basepath']
#     provider = references['provider']
#     release = references['release']

    # return [os.path.join(basepath, provider, release)]

def references_abs_path(ref='references'):
    references = config.get(ref)
    basepath = expand_filepath(references['basepath'])
    provider = references['provider']
    release = references['release']

    return [os.path.join(basepath, provider, release)]



def resolve_single_filepath(basepath, filename):
    return [os.path.join(basepath, filename)]

def resolve_multi_filepath(basepath, dictionary):
    for k, v in dictionary.items():
        dictionary[k] = os.path.join(basepath, v)
    return dictionary

def get_references_label(ref='references'):
    references = config.get(ref)
    provider = references['provider']
    genome = references['release']

    return '_'.join([provider, genome])





def get_units_by_sample(wildcards, samples, label='units', prefix='before',
                        suffix='after'):
    return [prefix+i+suffix for i in samples.loc[wildcards.sample,
                                                  [label]][0].split(',')]
def get_odp(wildcards,samples,optical_dup='odp'):
    return "OPTICAL_DUPLICATE_PIXEL_DISTANCE={}".format(samples.loc[wildcards.sample, [optical_dup]].dropna()[0])
###############################################################################

def get_sample_by_famid(wildcards, patients, label='sample'):
    for famid in patients.loc[wildcards.patient,[label]]:
        if samples.loc[famid.split(',')[0],"condition"]=="C" and samples.loc[famid.split(',')[1],"condition"]=="T":
#            print("all_ok")
            tbam = "reads/recalibrated/" + famid.split(',')[1] + ".dedup.recal.bam"
            cbam = "reads/recalibrated/" + famid.split(',')[0] + ".dedup.recal.bam"
        else:
            tbam = "reads/recalibrated/" + famid.split(',')[0] + ".dedup.recal.bam"
            cbam = "reads/recalibrated/" + famid.split(',')[1] + ".dedup.recal.bam"
    return (tbam,cbam)


def ensure_dir(path, force=False):
    try:
        if force and os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

