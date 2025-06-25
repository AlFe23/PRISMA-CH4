#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys

def load_image(tarball):
    print(f"▶️ Loading Docker image from {tarball} …")
    subprocess.run(['docker', 'load', '-i', tarball], check=True)

def ensure_dir(path):
    """If path is a file parent, return its dirname; 
       if it's a directory, ensure it exists and return it."""
    abs_path = os.path.abspath(path)
    if os.path.isdir(abs_path):
        os.makedirs(abs_path, exist_ok=True)
        return abs_path
    else:
        parent = os.path.dirname(abs_path)
        if not os.path.isdir(parent):
            os.makedirs(parent, exist_ok=True)
        return parent

def collect_mounts(paths):
    """Given a list of file/dir paths, return a dict {host_dir:host_dir} for mounting."""
    mounts = {}
    for p in paths:
        d = ensure_dir(p)
        mounts[d] = d
    return mounts

def build_docker_run_cmd(image, mounts, mode, script_args):
    cmd = ['docker', 'run', '--rm']
    for host_dir, cont_dir in mounts.items():
        cmd += ['-v', f'{host_dir}:{cont_dir}']
    cmd += [image, mode] + script_args
    return cmd

def main():
    parser = argparse.ArgumentParser(
        description="Load a CTMF Docker image from tarball and run it in single or batch mode"
    )
    parser.add_argument(
        '-t','--tarball',
        required=True,
        help="Path to your saved image tarball (e.g. prisma-ch4-mapper.tar)"
    )
    parser.add_argument(
        '--image',
        default='prisma-ch4-mapper:latest',
        help="Name:tag of the image inside the tarball (default: prisma-ch4-mapper:latest)"
    )

    sub = parser.add_subparsers(dest='mode', required=True, help='Choose single or batch')

    # single mode
    p1 = sub.add_parser('single', help='Process one L1+L2C acquisition')
    p1.add_argument('l1',    help='L1 .he5 file')
    p1.add_argument('l2c',   help='L2C .he5 file')
    p1.add_argument('dem',   help='DEM netCDF (.nc)')
    p1.add_argument('lut',   help='LUT HDF5 (.hdf5)')
    p1.add_argument('output',help='Output directory')
    p1.add_argument('-k', type=int, default=1, help='Cluster count')
    p1.add_argument('--min_wl', type=int, help='Min wavelength')
    p1.add_argument('--max_wl', type=int, help='Max wavelength')

    # batch mode
    p2 = sub.add_parser('batch', help='Process a whole directory of images')
    p2.add_argument('img_dir',     help='Directory of .he5 acquisitions')
    p2.add_argument('dem',         help='DEM netCDF (.nc)')
    p2.add_argument('lut',         help='LUT HDF5 (.hdf5)')
    p2.add_argument('--output_root', required=True, help='Output root directory')
    p2.add_argument('-k', type=int, default=1, help='Cluster count')
    p2.add_argument('--min_wl', type=int, help='Min wavelength')
    p2.add_argument('--max_wl', type=int, help='Max wavelength')

    args = parser.parse_args()

    # 1) load the image
    load_image(args.tarball)

    # 2) decide which paths to mount & which script-args to pass
    if args.mode == 'single':
        mounts = [args.l1, args.l2c, args.dem, args.lut, args.output]
        script_args = [
            args.l1, args.l2c, args.dem, args.lut, args.output,
            '-k', str(args.k)
        ]
        if args.min_wl: script_args += ['--min_wl', str(args.min_wl)]
        if args.max_wl: script_args += ['--max_wl', str(args.max_wl)]
    else:  # batch
        mounts = [args.img_dir, args.dem, args.lut, args.output_root]
        script_args = [
            args.img_dir, args.dem, args.lut,
            '--output_root', args.output_root,
            '-k', str(args.k)
        ]
        if args.min_wl: script_args += ['--min_wl', str(args.min_wl)]
        if args.max_wl: script_args += ['--max_wl', str(args.max_wl)]

    mount_map = collect_mounts(mounts)

    # 3) build & print the docker run command
    cmd = build_docker_run_cmd(args.image, mount_map, args.mode, script_args)
    print("▶️ Running container:\n   " + " ".join(cmd))

    # 4) actually run it
    subprocess.run(cmd, check=True)

if __name__ == '__main__':
    main()
