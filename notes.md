                        Traceback (most recent call last):
  File "/dev/shm/cybershuttle/envs/18dccede/bin/qcfractal-compute-manager", line 8, in 
    sys.exit(main())
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/site-packages/qcfractalcompute/compute_manager_cli.py", line 50, in main
    manager = ComputeManager(manager_config)
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/site-packages/qcfractalcompute/compute_manager.py", line 151, in __init__
    self.app_manager = AppManager(self.manager_config)
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/site-packages/qcfractalcompute/apps/app_manager.py", line 113, in __init__
    qcengine_functions = discover_programs_conda(conda_env)
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/site-packages/qcfractalcompute/apps/app_manager.py", line 33, in discover_programs_conda
    result = subprocess.check_output(cmd, universal_newlines=True, cwd=tmpdir)
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/subprocess.py", line 421, in check_output
    return run(*popenargs, stdout=PIPE, timeout=timeout, check=True,
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/subprocess.py", line 503, in run
    with Popen(*popenargs, **kwargs) as process:
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/subprocess.py", line 971, in __init__
    self._execute_child(args, executable, preexec_fn, close_fds,
  File "/dev/shm/cybershuttle/envs/18dccede/lib/python3.10/subprocess.py", line 1863, in _execute_child
    raise child_exception_type(errno_num, err_msg, err_filename)
FileNotFoundError: [Errno 2] No such file or directory: 'conda'
