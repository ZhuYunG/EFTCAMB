#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
import shutil
import re
from pathlib import Path
import numpy as np


# ===== 路径与基础设置 =====

# CAMB/EFTCosmoMC 可执行所在目录（按实际情况改）
CAMB_DIR = Path("/Users/dcz/Documents/eft_code/EFTCosmoMC/EFTCAMB")

# 模板 ini（完整 EFTCAMB 设置，Horndeski 参数会被覆盖）
TEMPLATE_INI = CAMB_DIR / "params_EFT.ini"

# 每次运行实际使用的 ini（脚本生成/覆盖）
RUN_INI = CAMB_DIR / "params.ini"

# 每次写入后的 ini 归档目录（可选）
RUN_INI_ARCHIVE_DIR = CAMB_DIR / "Horndeski_run_inis_py"

# Horndeski 样本输出目录
SAMPLES_DIR = CAMB_DIR / "Horndeski_samples_py"

# ./camb 可执行文件名（如果在 CAMB_DIR 下）
CAMB_EXE = CAMB_DIR / "camb"

# 采样次数
N_SAMPLES = 10000


# ===== 1. 设定 Horndeski 参数及其采样范围 =====
# 这里用一个字典描述参数名及其采样区间，名字必须和 ini 文件中的键一致
PARAM_RANGES = {
    #根据你自己的 Horndeski 模块的 parameter_names 来填
    "OmegaN0": (-1.0, 1.0),
    "OmegaN1": (-1.0, 1.0),
    "OmegaN2": (-1.0, 1.0),
    "OmegaN3": (-1.0, 1.0),
    "OmegaN4": (-1.0, 1.0),
    "OmegaN5": (-1.0, 1.0),
    "OmegaD1": (-1.0, 1.0),
    "OmegaD2": (-1.0, 1.0),
    "OmegaD3": (-1.0, 1.0),
    "OmegaD4": (-1.0, 1.0),
    "LambdaN0": (-1.0, 1.0),
    "LambdaN1": (-1.0, 1.0),
    "LambdaN2": (-1.0, 1.0),
    "LambdaN3": (-1.0, 1.0),
    "LambdaN4": (-1.0, 1.0),
    "LambdaN5": (-1.0, 1.0),
    "LambdaD1": (-1.0, 1.0),
    "LambdaD2": (-1.0, 1.0),
    "LambdaD3": (-1.0, 1.0),
    "LambdaD4": (-1.0, 1.0)
}


def sample_parameters():
    """在 PARAM_RANGES 规定的区间内均匀采样一组参数。"""
    params = {}
    for name, (pmin, pmax) in PARAM_RANGES.items():
        u = np.random.rand()
        params[name] = pmin + u * (pmax - pmin)
    return params


# ===== 2. 从模板 ini 生成本次运行用的 ini =====

def update_ini_text(template_text, new_params, file_root):
    """
    用 new_params 覆盖模板中的 Horndeski 参数行，
    同时设置 file_root = file_root。
    CosmoMC/EFTCAMB 的 ini 格式一般是 'name = value' 或 'name=value'。
    """
    lines = template_text.splitlines()
    out_lines = []

    # 方便匹配：键名 -> (pattern, new_line)
    patterns = {}
    for name, value in new_params.items():
        # 正则：行首可有空格，然后是 name，然后是可选空格、= 号
        pat = re.compile(rf"^\s*{re.escape(name)}\s*=")
        # 统一写成 name = value 格式
        new_line = f"{name} = {value:.8e}"
        patterns[name] = (pat, new_line)

    # file_root 也单独处理
    pat_file_root = re.compile(r"^\s*file_root\s*=")
    new_file_root_line = f"file_root = {file_root}"

    replaced_params = set()
    file_root_replaced = False

    for line in lines:
        stripped = line.strip()

        # 注释或空行不动
        if not stripped or stripped.startswith("#"):
            out_lines.append(line)
            continue

        # 先处理 file_root
        if pat_file_root.match(line):
            out_lines.append(new_file_root_line)
            file_root_replaced = True
            continue

        # 再处理所有 Horndeski 参数
        replaced = False
        for name, (pat, new_line) in patterns.items():
            if pat.match(line):
                out_lines.append(new_line)
                replaced = True
                replaced_params.add(name)
                break

        if not replaced:
            out_lines.append(line)

    # 将未出现过的参数追加到文件末尾，确保每个参数都有覆盖
    for name, (_, new_line) in patterns.items():
        if name not in replaced_params:
            out_lines.append(new_line)

    # file_root 如果模板中不存在，则追加
    if not file_root_replaced:
        out_lines.append(new_file_root_line)

    return "\n".join(out_lines) + "\n"


def write_run_ini(params, sample_id):
    """
    从 params_EFT.ini 生成本次样本的配置，覆盖写回 params_EFT.ini，
    并归档一份；params.ini 不做改动（camb 仍使用现有的 params.ini）。
    """
    if not TEMPLATE_INI.is_file():
        raise FileNotFoundError(f"找不到模板 ini: {TEMPLATE_INI}")

    text = TEMPLATE_INI.read_text(encoding="utf-8")

    # 为每个样本生成一个独立的 file_root，避免输出覆盖
    file_root = f"hsamp_{sample_id:05d}"

    new_text = update_ini_text(text, params, file_root)

    # 覆盖写 params_EFT.ini（用户指定）
    TEMPLATE_INI.write_text(new_text, encoding="utf-8")

    # 归档本次 ini，方便回溯
    RUN_INI_ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    archive_path = RUN_INI_ARCHIVE_DIR / f"params_EFT_{sample_id:05d}.ini"
    archive_path.write_text(new_text, encoding="utf-8")

    return file_root


# ===== 3. 调用 ./camb 运行一次 EFTCAMB =====

def run_camb(run_ini_path):
    """
    调用 ./camb run_ini_path。
    返回 (success, stdout, stderr)。
    success 仅根据 returncode 判断，进一步的物理稳定性可以在 stdout 里再查。
    """
    cmd = [str(CAMB_EXE), str(run_ini_path)]
    # 在 CAMB_DIR 下运行，保持相对路径一致
    result = subprocess.run(
        cmd,
        cwd=str(CAMB_DIR),
        capture_output=True,
        text=True
    )

    success = (result.returncode == 0)
    return success, result.stdout, result.stderr


# ===== 4. 根据 file_root 找到背景输出并保存为 Horndeski_sample_i.dat =====

# def collect_background(file_root, sample_id):
#     """
#     根据 file_root 找到背景输出文件，复制/重命名到 Horndeski_samples/Horndeski_sample_{sample_id}.dat。
#     具体文件名要根据你 EFTCAMB 的实际输出调整。
#     """
#     SAMPLES_DIR.mkdir(parents=True, exist_ok=True)

#     # 这里是假设背景文件名形如  file_root + "_background.dat"
#     # 你需要根据自己实际的输出文件名来改，比如 "_eftback.dat" 等。
#     background_name = f"{file_root}_background.dat"
#     background_path = CAMB_DIR / background_name

#     if not background_path.is_file():
#         raise FileNotFoundError(f"找不到背景输出文件: {background_path}")

#     target_path = SAMPLES_DIR / f"Horndeski_sample_{sample_id}.dat"
#     shutil.copy2(background_path, target_path)
#     return target_path

def collect_background(sample_id):
    """
    从固定路径 Horndeski_solution.dat 复制为 Horndeski_sample_{sample_id}.dat
    """
    SAMPLES_DIR.mkdir(parents=True, exist_ok=True)

    # 固定输出文件
    background_path = CAMB_DIR / "Horndeski_solution.dat"

    if not background_path.is_file():
        raise FileNotFoundError(f"找不到背景输出文件: {background_path}")

    target_path = SAMPLES_DIR / f"Horndeski_sample_{sample_id}.dat"
    shutil.copy2(background_path, target_path)
    return target_path



# ===== 5. 主循环：采样 → 写 ini → 跑 camb → 检查 → 收集 =====

def main():
    if not PARAM_RANGES:
        raise RuntimeError("PARAM_RANGES 为空，请先在脚本顶部填入 Horndeski 参数名及采样范围。")

    if not CAMB_EXE.is_file():
        raise FileNotFoundError(f"找不到可执行文件: {CAMB_EXE}")

    print(f"Working dir: {CAMB_DIR}")
    print(f"Template ini: {TEMPLATE_INI}")
    print(f"Samples dir: {SAMPLES_DIR}")
    print(f"Total samples: {N_SAMPLES}")

    success_count = 0

    for i in range(1, N_SAMPLES + 1):
        print(f"\n=== Sample {i}/{N_SAMPLES} ===")

        # 1) 采样参数
        params = sample_parameters()
        print("  parameters:", params)

        # 2) 写 params_EFT_run.ini 并设置 file_root
        file_root = write_run_ini(params, i)
        print("  file_root:", file_root)

        # 3) 调用 ./camb
        success, out, err = run_camb(RUN_INI)

        # 4) 用 returncode 做第一次筛选
        if not success:
            print("  CAMB run failed, returncode != 0")
            print("  stderr (first lines):")
            for line in err.splitlines()[:5]:
                print("    ", line)
            continue

        print("  CAMB stdout (first lines):")
        for line in out.splitlines()[:5]:
            print("    ", line)

        # 5) 可选：在 stdout 中搜索不稳定关键字，再进一步筛掉
        # 例如：
        # if "EFTCAMB ghost instability" in out or "EFTCAMB gradient instability" in out:
        #     print("  EFTCAMB signaled instability, reject sample.")
        #     continue

        # 6) 收集背景文件
        # try:
        #     target_path = collect_background(file_root, i)
        #     print(f"  background saved to {target_path}")
        #     success_count += 1
        # except FileNotFoundError as e:
        #     print("  WARNING:", e)
        #     # 如果背景文件不存在，可以把这次样本也视为失败
        #     continue

        # 6) 收集背景文件（从 Horndeski_solution.dat 复制）
        try:
            target_path = collect_background(i)
            print(f"  background saved to {target_path}")
            success_count += 1
        except FileNotFoundError as e:
            print("  WARNING:", e)
            continue


    print(f"\nDone. Accepted {success_count} / {N_SAMPLES} samples.")


if __name__ == "__main__":
    main()
