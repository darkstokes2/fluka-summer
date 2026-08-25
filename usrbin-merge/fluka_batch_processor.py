#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# fluka-summer — FLUKA post-processing toolkit
# Copyright (C) 2026  Pengbin Zhang
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
批处理脚本：用 usbsuw 合并各文件夹的 _fort.42 文件，再用 usbrea 转成 .lis
- 输入: 各文件夹下的 *_fort.NN (命名规律: inpX00Y_fort.NN)
- usbsuw: 多个 _fort.NN → 1个 .bnn (放各自文件夹)
- usbrea: .bnn → .lis (放 postprocess 文件夹)

用法:
  # 自动发现 base 目录下含 _fort.42 的子文件夹
  python fluka_batch_processor.py --base /path/to/simulation/projects

  # 指定文件夹列表和 fort 单元号
  python fluka_batch_processor.py --base /path/to/projects \\
      --folders run_A run_B --fort-unit 42 --tag Proton

  # 指定 FLUKA 工具路径和超时
  python fluka_batch_processor.py --base /path \\
      --usbsuw-bin /opt/fluka/bin/usbsuw \\
      --usbrea-bin /opt/fluka/bin/usbrea \\
      --usbsuw-timeout 600 --usbrea-timeout 1200
"""

import argparse
import glob
import os
import sys
import time

import pexpect


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch-process FLUKA USRBIN outputs: "
                    "merge _fort.NN files with usbsuw, convert to .lis with usbrea.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--base", required=True,
        help="Base directory containing simulation sub-folders.",
    )
    parser.add_argument(
        "--post", default=None,
        help="Output directory for .lis files (default: <base>/postprocess).",
    )
    parser.add_argument(
        "--folders", nargs="*", default=None,
        help="Sub-folder names to process (default: auto-discover folders "
             "containing _fort.<unit> files).",
    )
    parser.add_argument(
        "--fort-unit", default="42",
        help="Fortran unit number in _fort.NN (default: 42 = energy deposition).",
    )
    parser.add_argument(
        "--tag", default="Proton",
        help="Label for output file naming (default: Proton).",
    )
    parser.add_argument(
        "--usbsuw-bin", default="usbsuw",
        help="Path to usbsuw executable (default: usbsuw from PATH).",
    )
    parser.add_argument(
        "--usbrea-bin", default="usbrea",
        help="Path to usbrea executable (default: usbrea from PATH).",
    )
    parser.add_argument(
        "--usbsuw-timeout", type=int, default=300,
        help="Timeout in seconds for usbsuw merge (default: 300).",
    )
    parser.add_argument(
        "--usbrea-timeout", type=int, default=600,
        help="Timeout in seconds for usbrea conversion (default: 600).",
    )
    return parser.parse_args()


def discover_folders(base, fort_unit):
    """Auto-discover sub-folders containing _fort.<unit> files."""
    folders = []
    for entry in sorted(os.listdir(base)):
        folder_path = os.path.join(base, entry)
        if not os.path.isdir(folder_path):
            continue
        pattern = os.path.join(folder_path, f"*_fort.{fort_unit}")
        if glob.glob(pattern):
            folders.append(entry)
    if not folders:
        print(f"[Warning] No folders containing _fort.{fort_unit} found in {base}")
    return folders


def run_usbsuw(usbsuw_bin, folder_path, fort_files, bnn_name, timeout):
    """用 pexpect 驱动 usbsuw，合并多个 fort 文件为一个 bnn"""
    child = pexpect.spawn(usbsuw_bin, encoding="utf-8", timeout=timeout)
    child.logfile_read = sys.stdout  # 实时打印输出

    for fort in fort_files:
        child.expect("Type the input file:")
        child.sendline(fort)

    # 空回车结束输入
    child.expect("Type the input file:")
    child.sendline("")

    # 输出文件名
    child.expect("Type the output file name:")
    child.sendline(bnn_name)

    child.expect(pexpect.EOF)
    child.close()
    return os.path.exists(os.path.join(folder_path, bnn_name))


def run_usbrea(usbrea_bin, folder_path, bnn_name, lis_path, timeout):
    """用 pexpect 驱动 usbrea，把 bnn 转成 lis"""
    child = pexpect.spawn(usbrea_bin, encoding="utf-8", timeout=timeout)
    child.logfile_read = sys.stdout

    bnn_path = os.path.join(folder_path, bnn_name)
    child.expect("Type the input file:")
    child.sendline(bnn_path)

    child.expect("Type the output file name:")
    child.sendline(lis_path)

    child.expect(pexpect.EOF)
    child.close()
    return os.path.exists(lis_path)


def process_folder(folder, base, post, fort_unit, tag,
                   usbsuw_bin, usbrea_bin, usbsuw_timeout, usbrea_timeout):
    """处理单个文件夹"""
    folder_path = os.path.join(base, folder)
    pattern = os.path.join(folder_path, f"*_fort.{fort_unit}")
    fort_files = sorted(glob.glob(pattern))

    if len(fort_files) == 0:
        print(f"\n[{folder}] 无 _fort.{fort_unit} 文件，跳过")
        return False

    print(f"\n{'='*60}")
    print(f"[{folder}] 找到 {len(fort_files)} 个 _fort.{fort_unit} 文件")
    print(f"  第一个: {os.path.basename(fort_files[0])}")
    print(f"  最后一个: {os.path.basename(fort_files[-1])}")

    # 文件名用相对路径（cd 到文件夹后执行）
    fort_rel = [os.path.basename(f) for f in fort_files]
    bnn_name = f"{folder}_{tag}.bnn"
    bnn_path = os.path.join(folder_path, bnn_name)
    lis_name = f"{folder}_{tag}.lis"
    lis_path = os.path.join(post, lis_name)

    # 步骤1: usbsuw 合并
    print(f"\n--- 步骤1: usbsuw 合并 {len(fort_rel)} 个文件 → {bnn_name} ---")
    t0 = time.time()
    # usbsuw 需要在目标文件夹下运行（文件名带空格/括号，用相对路径最安全）
    os.chdir(folder_path)
    ok = run_usbsuw(usbsuw_bin, folder_path, fort_rel, bnn_name, usbsuw_timeout)
    t1 = time.time()
    if not ok:
        print(f"  [失败] usbsuw 未生成 {bnn_name}")
        return False
    size_mb = os.path.getsize(bnn_path) / 1e6
    print(f"  [完成] {bnn_name} ({size_mb:.1f} MB, 用时 {t1-t0:.1f}s)")

    # 步骤2: usbrea 转换
    print(f"\n--- 步骤2: usbrea 转换 → {lis_name} ---")
    t0 = time.time()
    ok = run_usbrea(usbrea_bin, folder_path, bnn_name, lis_path, usbrea_timeout)
    t1 = time.time()
    if not ok:
        print(f"  [失败] usbrea 未生成 {lis_name}")
        return False
    size_mb = os.path.getsize(lis_path) / 1e6
    print(f"  [完成] {lis_name} ({size_mb:.1f} MB, 用时 {t1-t0:.1f}s)")

    return True


def main():
    args = parse_args()

    base = args.base
    post = args.post or os.path.join(base, "postprocess")
    fort_unit = args.fort_unit
    tag = args.tag

    # 文件夹列表：优先用命令行参数，否则自动发现
    if args.folders:
        folders = args.folders
    else:
        folders = discover_folders(base, fort_unit)

    os.makedirs(post, exist_ok=True)
    print(f"基础路径: {base}")
    print(f"输出路径: {post}")
    print(f"处理单元: _fort.{fort_unit} ({tag})")
    print(f"文件夹数: {len(folders)}")
    print(f"usbsuw: {args.usbsuw_bin}")
    print(f"usbrea: {args.usbrea_bin}")

    original_cwd = os.getcwd()
    results = {}
    for folder in folders:
        ok = process_folder(
            folder, base, post, fort_unit, tag,
            args.usbsuw_bin, args.usbrea_bin,
            args.usbsuw_timeout, args.usbrea_timeout,
        )
        results[folder] = ok
    os.chdir(original_cwd)

    # 汇总
    print(f"\n{'='*60}")
    print("=== 最终汇总 ===")
    for folder, ok in results.items():
        status = "✅ 成功" if ok else "❌ 失败"
        print(f"  {folder:15s} {status}")

    # 列出生成的 .lis 文件
    print(f"\n=== postprocess 文件夹下的 .lis 文件 ===")
    for f in sorted(glob.glob(os.path.join(post, "*.lis"))):
        size_mb = os.path.getsize(f) / 1e6
        print(f"  {os.path.basename(f):40s} {size_mb:.1f} MB")


if __name__ == "__main__":
    main()
