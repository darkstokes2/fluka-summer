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
批处理脚本：用 usxsuw 合并各文件夹的 USRBDX _fort.NN 文件
- 输入: 各文件夹下的 *_fort.NN (命名规律: inpX00Y_fort.NN)
- usxsuw: 多个 _fort.NN → 1个 .bnn + _sum.lis + _tab.lis
- .bnn 放各自文件夹, _sum.lis 和 _tab.lis 放 postprocess 文件夹
- usxsuw 会自动生成3个文件: xxx.bnn / xxx_sum.lis / xxx_tab.lis (无需额外调 usbrea)

用法:
  # 自动发现 base 目录下含 _fort.50 的子文件夹，处理 Proton 探测器
  python usrbdx_batch_processor.py --base /path/to/simulation/projects

  # 指定文件夹列表，处理多个探测器
  python usrbdx_batch_processor.py --base /path/to/projects \\
      --folders run_A run_B \\
      --fort-units 50 51 52 --tags Proton Neutron Electron

  # 指定 FLUKA 工具路径和超时
  python usrbdx_batch_processor.py --base /path \\
      --usxsuw-bin /opt/fluka/bin/usxsuw --usxsuw-timeout 600
"""

import argparse
import glob
import os
import sys
import time

import pexpect


def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch-process FLUKA USRBDX outputs: "
                    "merge _fort.NN files with usxsuw into .bnn, _sum.lis, _tab.lis.",
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
        "--fort-units", nargs="+", default=["50"],
        help="USRBDX Fortran unit numbers (default: 50). "
             "Example: --fort-units 50 51 52",
    )
    parser.add_argument(
        "--tags", nargs="+", default=["Proton"],
        help="Labels for each detector, parallel to --fort-units (default: Proton). "
             "Example: --tags Proton Neutron Electron",
    )
    parser.add_argument(
        "--usxsuw-bin", default="usxsuw",
        help="Path to usxsuw executable (default: usxsuw from PATH).",
    )
    parser.add_argument(
        "--usxsuw-timeout", type=int, default=300,
        help="Timeout in seconds per usxsuw merge (default: 300).",
    )
    return parser.parse_args()


def discover_folders(base, fort_units):
    """Auto-discover sub-folders containing _fort.<unit> files for any unit."""
    folders = []
    for entry in sorted(os.listdir(base)):
        folder_path = os.path.join(base, entry)
        if not os.path.isdir(folder_path):
            continue
        for unit in fort_units:
            pattern = os.path.join(folder_path, f"*_fort.{unit}")
            if glob.glob(pattern):
                folders.append(entry)
                break
    if not folders:
        print(f"[Warning] No folders containing _fort.{'/'.join(fort_units)} "
              f"found in {base}")
    return folders


def run_usxsuw(usxsuw_bin, folder_path, fort_files, bnn_name, timeout):
    """用 pexpect 驱动 usxsuw，合并多个 fort 文件为一个 bnn (+ _sum.lis + _tab.lis)"""
    child = pexpect.spawn(usxsuw_bin, encoding="utf-8", timeout=timeout)
    child.logfile_read = sys.stdout

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
    # usxsuw 生成3个文件: xxx.bnn, xxx_sum.lis, xxx_tab.lis
    return os.path.exists(os.path.join(folder_path, bnn_name))


def process_unit(folder, fort_unit, tag, base, post,
                 usxsuw_bin, usxsuw_timeout):
    """处理单个文件夹的单个 fort 单元"""
    folder_path = os.path.join(base, folder)
    pattern = os.path.join(folder_path, f"*_fort.{fort_unit}")
    fort_files = sorted(glob.glob(pattern))

    if len(fort_files) == 0:
        print(f"  [{folder}/fort.{fort_unit}] 无文件，跳过")
        return False

    print(f"\n  [fort.{fort_unit} ({tag})] {len(fort_files)} 个文件")

    # 文件名用相对路径（cd 到文件夹后执行）
    fort_rel = [os.path.basename(f) for f in fort_files]
    bnn_name = f"{folder}_bdx_{tag}.bnn"

    # 步骤: usxsuw 合并（会自动生成 _sum.lis 和 _tab.lis）
    t0 = time.time()
    os.chdir(folder_path)
    ok = run_usxsuw(usxsuw_bin, folder_path, fort_rel, bnn_name, usxsuw_timeout)
    t1 = time.time()
    if not ok:
        print(f"    [失败] usxsuw 未生成 {bnn_name}")
        return False
    print(f"    [完成] {bnn_name} ({os.path.getsize(bnn_name)/1e3:.1f} KB, {t1-t0:.1f}s)")

    # 把 _sum.lis 和 _tab.lis 移到 postprocess
    base_no_ext = bnn_name[:-4]  # 去掉 .bnn
    for suffix in ["_sum.lis", "_tab.lis"]:
        src = os.path.join(folder_path, base_no_ext + suffix)
        dst = os.path.join(post, base_no_ext + suffix)
        if os.path.exists(src):
            # 移动而非复制（.bnn 留在文件夹里，.lis 移走）
            os.rename(src, dst)
            print(f"    → {os.path.basename(dst)} ({os.path.getsize(dst)/1e3:.1f} KB)")

    return True


def main():
    args = parse_args()

    base = args.base
    post = args.post or os.path.join(base, "postprocess")
    fort_units = args.fort_units
    tags = args.tags

    if len(fort_units) != len(tags):
        print("[Error] --fort-units 和 --tags 的数量必须一致")
        sys.exit(1)

    # 文件夹列表：优先用命令行参数，否则自动发现
    if args.folders:
        folders = args.folders
    else:
        folders = discover_folders(base, fort_units)

    os.makedirs(post, exist_ok=True)
    print(f"基础路径: {base}")
    print(f"输出路径: {post}")
    print(f"处理单元: {fort_units} (USRBDX)")
    print(f"文件夹数: {len(folders)}")
    print(f"usxsuw: {args.usxsuw_bin}")

    original_cwd = os.getcwd()
    results = {}
    for folder in folders:
        print(f"\n{'='*60}")
        print(f"[{folder}]")
        folder_results = []
        for fort_unit, tag in zip(fort_units, tags):
            ok = process_unit(folder, fort_unit, tag, base, post,
                              args.usxsuw_bin, args.usxsuw_timeout)
            folder_results.append(ok)
        results[folder] = folder_results
    os.chdir(original_cwd)

    # 汇总
    print(f"\n{'='*60}")
    print("=== 最终汇总 ===")
    for folder, oks in results.items():
        for tag, ok in zip(tags, oks):
            status = "✅" if ok else "❌"
            print(f"  {folder:15s} fort.{fort_units[tags.index(tag)]} {tag:10s} {status}")

    # 列出生成的文件
    print(f"\n=== postprocess 文件夹下的 USRBDX 文件 ===")
    for f in sorted(glob.glob(os.path.join(post, "*_bdx_*.lis"))):
        size_kb = os.path.getsize(f) / 1e3
        print(f"  {os.path.basename(f):45s} {size_kb:.1f} KB")


if __name__ == "__main__":
    main()
