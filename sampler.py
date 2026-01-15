#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import json
import os
import subprocess
import shutil
import re
import time
import threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path
import numpy as np


# ===== 路径与基础设置 =====

# CAMB/EFTCosmoMC 可执行所在目录（按实际情况改）
CAMB_DIR = Path("/Users/dcz/Documents/eft_code/EFTCosmoMC/EFTCAMB")

# 模板 ini（完整 EFTCAMB 设置，Horndeski 参数会被覆盖）
TEMPLATE_EFT_INI = CAMB_DIR / "params_EFT.ini"

# 运行 ini 模板（会在运行时复制生成每个样本的独立 ini）
TEMPLATE_RUN_INI = CAMB_DIR / "params.ini"

# 每次写入后的 ini 归档目录（可选，可通过外部参数/环境变量覆盖）
DEFAULT_RUN_INI_ARCHIVE_DIR = Path("/Users/dcz/data/Horndeski_samples/Horndeski_run_inis_py_onlybackground_a01_2")

# Horndeski 样本输出目录（可通过外部参数/环境变量覆盖）
DEFAULT_SAMPLES_DIR = Path("/Users/dcz/data/Horndeski_samples/Horndeski_samples_py_onlybackground_a01_2")

# 每个样本运行的独立工作目录基路径（可通过外部参数/环境变量覆盖）
DEFAULT_RUN_WORK_DIR_BASE = Path("/Volumes/My Passport/Horndeski_workdirs_py_1")

# ./camb 可执行文件名（如果在 CAMB_DIR 下）
CAMB_EXE = CAMB_DIR / "camb"

# 采样次数
N_SAMPLES = 250000

# 单个样本 CAMB 运行超时（秒），超过则判为失败
CAMB_TIMEOUT_SECONDS = 6 * 60


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
    "LambdaD4": (-1.0, 1.0),
    "hubble": (20.0, 100.0),
    "omegaLambda": (0.0, 1.0)
}

BASE_OMBH2 = 0.0226
BASE_OMCH2 = 0.112
F_B = BASE_OMBH2 / (BASE_OMBH2 + BASE_OMCH2)   # ~0.1679049

OMEGA_M_RANGE = (0.0, 1.0)  # Omega_m in [0,1]


def parse_args():
    parser = argparse.ArgumentParser(description="EFTCAMB sampler")
    parser.add_argument(
        "--run-ini-archive-dir",
        default=None,
        help="Archive dir for each generated ini (or set RUN_INI_ARCHIVE_DIR)."
    )
    parser.add_argument(
        "--save-run-ini",
        dest="save_run_ini",
        action="store_true",
        help="Keep and archive per-sample ini files (or set SAVE_RUN_INI=1)."
    )
    parser.add_argument(
        "--no-save-run-ini",
        dest="save_run_ini",
        action="store_false",
        help="Do not keep/archive per-sample ini files."
    )
    parser.add_argument(
        "--samples-dir",
        default=None,
        help="Output dir for Horndeski samples (or set SAMPLES_DIR)."
    )
    parser.add_argument(
        "--run-work-dir-base",
        default=None,
        help="Base dir for per-sample working dirs (or set RUN_WORK_DIR_BASE)."
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=None,
        help="Parallel workers for CAMB runs (or set MAX_WORKERS)."
    )
    parser.add_argument(
        "--omp-threads",
        type=int,
        default=None,
        help="OMP_NUM_THREADS value for CAMB (defaults to 1 when max-workers > 1)."
    )
    parser.add_argument(
        "--sample-id-start",
        type=int,
        default=None,
        help="Start index for output sample filenames (or set SAMPLE_ID_START)."
    )
    parser.add_argument(
        "--parallel",
        dest="parallel",
        action="store_true",
        help="Enable parallel CAMB runs (or set USE_PARALLEL=1)."
    )
    parser.add_argument(
        "--no-parallel",
        dest="parallel",
        action="store_false",
        help="Disable parallel mode and use legacy sequential behavior."
    )
    parser.set_defaults(parallel=None, save_run_ini=None)
    return parser.parse_args()


def resolve_dir(cli_value, env_key, default_path):
    if cli_value:
        return Path(cli_value).expanduser()
    env_value = os.environ.get(env_key)
    if env_value:
        return Path(env_value).expanduser()
    return default_path


def resolve_int(cli_value, env_key, default_value):
    if cli_value is not None:
        return int(cli_value)
    env_value = os.environ.get(env_key)
    if env_value:
        return int(env_value)
    return default_value


def resolve_bool(cli_value, env_key, default_value):
    if cli_value is not None:
        return bool(cli_value)
    env_value = os.environ.get(env_key)
    if env_value is None:
        return default_value
    value = env_value.strip().lower()
    if value in {"1", "true", "yes", "y", "on"}:
        return True
    if value in {"0", "false", "no", "n", "off"}:
        return False
    raise ValueError(f"{env_key} 取值无效: {env_value}")


def sample_parameters(rng=None):
    """在 PARAM_RANGES 规定的区间内均匀采样一组参数。"""
    rng = rng or np.random
    params = {}
    for name, (pmin, pmax) in PARAM_RANGES.items():
        u = rng.random()
        params[name] = pmin + u * (pmax - pmin)

    Omega_m = rng.uniform(*OMEGA_M_RANGE)

    # 用 hubble 计算 h，并把 Omega_m -> (ombh2, omch2)
    H0 = params["hubble"]          # 你的 PARAM_RANGES 里已有 hubble
    h = H0 / 100.0
    omega_m = Omega_m * (h * h)    # omega_m = Omega_m * h^2

    params["ombh2"] = F_B * omega_m
    params["omch2"] = (1.0 - F_B) * omega_m

    return params, Omega_m


def init_sample_log(samples_dir):
    samples_dir.mkdir(parents=True, exist_ok=True)
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    log_path = samples_dir / f"sample_log_{timestamp}_{time.time_ns()}.txt"
    header = (
        "sample_id\tsuccess\thubble\tomegaLambda\tomega_m\tOmega_m_draw\tparams_json\n"
    )
    log_path.write_text(header, encoding="utf-8")
    return log_path


def format_range_value(value):
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (int, float)):
        return repr(value)
    return str(value)


def format_range(range_value):
    if isinstance(range_value, (tuple, list)) and len(range_value) == 2:
        return f"[{format_range_value(range_value[0])}, {format_range_value(range_value[1])}]"
    return str(range_value)


def append_run_readme(samples_dir, sample_id_start):
    samples_dir.mkdir(parents=True, exist_ok=True)
    run_id = f"{time.strftime('%Y%m%d_%H%M%S')}_{time.time_ns()}"
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    lines = [
        f"=== Run {run_id} ===",
        f"Timestamp: {timestamp}",
        f"N_SAMPLES: {N_SAMPLES}",
        f"SAMPLE_ID_START: {sample_id_start}",
        f"BASE_OMBH2: {format_range_value(BASE_OMBH2)}",
        f"BASE_OMCH2: {format_range_value(BASE_OMCH2)}",
        f"OMEGA_M_RANGE: {format_range(OMEGA_M_RANGE)}",
        "",
        "PARAM_RANGES:",
    ]
    for name, rng in PARAM_RANGES.items():
        lines.append(f"  {name}: {format_range(rng)}")
    lines.append("Progress: 0 / 0 (running)")
    readme_path = samples_dir / "README.txt"
    separator = ""
    if readme_path.is_file() and readme_path.stat().st_size > 0:
        separator = "\n"
    with readme_path.open("a", encoding="utf-8") as handle:
        if separator:
            handle.write(separator)
        handle.write("\n".join(lines) + "\n")
    return readme_path, run_id


def update_run_readme_progress(samples_dir, run_id, accepted, processed, status):
    readme_path = samples_dir / "README.txt"
    if not readme_path.is_file():
        return None
    lines = readme_path.read_text(encoding="utf-8").splitlines()
    header = f"=== Run {run_id} ==="
    header_index = None
    for i, line in enumerate(lines):
        if line.strip() == header:
            header_index = i
            break
    if header_index is None:
        return None
    end_index = len(lines)
    for i in range(header_index + 1, len(lines)):
        if lines[i].startswith("=== Run "):
            end_index = i
            break
    progress_line = f"Progress: {accepted} / {processed} ({status})"
    updated = False
    for i in range(header_index + 1, end_index):
        if lines[i].startswith("Progress:"):
            lines[i] = progress_line
            updated = True
            break
    if not updated:
        lines.insert(end_index, progress_line)
    readme_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return readme_path


def compute_omega_m(params, omega_m_draw):
    omega_m = None
    if "ombh2" in params and "omch2" in params:
        try:
            omega_m = float(params["ombh2"]) + float(params["omch2"])
        except (TypeError, ValueError):
            omega_m = None
    if omega_m is None:
        hubble = params.get("hubble")
        if hubble is not None and omega_m_draw is not None:
            try:
                h = float(hubble) / 100.0
                omega_m = float(omega_m_draw) * (h * h)
            except (TypeError, ValueError):
                omega_m = None
    return omega_m


def normalize_params(params):
    normalized = {}
    for key, value in params.items():
        if isinstance(value, np.generic):
            normalized[key] = value.item()
        elif isinstance(value, (int, float)):
            normalized[key] = float(value)
        else:
            normalized[key] = str(value)
    return normalized


def format_field(value):
    if value is None:
        return ""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, (int, float)):
        return f"{float(value):.8e}"
    return str(value)


def append_sample_log(log_path, log_lock, sample_id, params, omega_m_draw, success):
    hubble = params.get("hubble")
    omega_lambda = params.get("omegaLambda")
    omega_m = compute_omega_m(params, omega_m_draw)
    params_json = json.dumps(
        normalize_params(params),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    fields = [
        str(sample_id),
        "1" if success else "0",
        format_field(hubble),
        format_field(omega_lambda),
        format_field(omega_m),
        format_field(omega_m_draw),
        params_json,
    ]
    line = "\t".join(fields) + "\n"
    if log_lock is not None:
        with log_lock:
            with log_path.open("a", encoding="utf-8") as handle:
                handle.write(line)
    else:
        with log_path.open("a", encoding="utf-8") as handle:
            handle.write(line)

# 以下用于调试排查pipeline本身是否存在问题, 根据测试, pipeline本身没有问题
# def sample_parameters():
#     # 暂时用一组你知道能跑通的参数，先排查管线
#     return {
#         "OmegaN0": 0.00,
#         "OmegaN1": 0.00,
#         "OmegaN2": 0.00,
#         "OmegaN3": 0.00,
#         "OmegaN4": 0.00,
#         "OmegaN5": 0.00,
#         "OmegaD1": 0.00,
#         "OmegaD2": 0.00,
#         "OmegaD3": 0.00,
#         "OmegaD4": 0.00,
#         "LambdaN0": 0.00,
#         "LambdaN1": 0.00,
#         "LambdaN2": 0.00,
#         "LambdaN3": 0.00,
#         "LambdaN4": 0.00,
#         "LambdaN5": 0.00,
#         "LambdaD1": 0.00,
#         "LambdaD2": 0.00,
#         "LambdaD3": 0.00,
#         "LambdaD4": 0.00,
#     }



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


def update_run_ini_text(template_text, eft_ini_name, highl_template_path):
    """
    生成 params.ini 的运行版本，确保 DEFAULT(...) 指向本次样本的 EFT ini，
    并固定 highL_unlensed_cl_template 为绝对路径，避免工作目录变更导致找不到模板。
    """
    lines = template_text.splitlines()
    out_lines = []

    pat_default = re.compile(r"^\s*DEFAULT\s*\(")
    pat_highl = re.compile(r"^\s*highL_unlensed_cl_template\s*=")

    default_replaced = False
    highl_replaced = False

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            out_lines.append(line)
            continue

        if pat_default.match(line):
            out_lines.append(f"DEFAULT({eft_ini_name})")
            default_replaced = True
            continue

        if pat_highl.match(line):
            out_lines.append(f"highL_unlensed_cl_template = {highl_template_path}")
            highl_replaced = True
            continue

        out_lines.append(line)

    if not default_replaced:
        out_lines.insert(0, f"DEFAULT({eft_ini_name})")

    if not highl_replaced:
        out_lines.append(f"highL_unlensed_cl_template = {highl_template_path}")

    return "\n".join(out_lines) + "\n"


def make_run_dir(run_work_base_dir, sample_id):
    run_dir = run_work_base_dir / f"run_{sample_id:05d}"
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        run_dir = run_work_base_dir / f"run_{sample_id:05d}_{time.time_ns()}"
        run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def write_run_files(
    params,
    sample_id,
    run_dir,
    run_ini_archive_dir,
    save_run_ini,
    eft_template_text,
    run_template_text,
    highl_template_path,
):
    """
    为每个样本生成独立的 params_EFT_XXXXX.ini 与 params_XXXXX.ini。
    当 save_run_ini=True 时归档 EFT ini。
    """
    file_root = f"hsamp_{sample_id:05d}"

    eft_text = update_ini_text(eft_template_text, params, file_root)
    eft_ini_path = run_dir / f"params_EFT_{sample_id:05d}.ini"
    eft_ini_path.write_text(eft_text, encoding="utf-8")

    if save_run_ini and run_ini_archive_dir is not None:
        run_ini_archive_dir.mkdir(parents=True, exist_ok=True)
        archive_path = run_ini_archive_dir / f"params_EFT_{sample_id:05d}.ini"
        archive_path.write_text(eft_text, encoding="utf-8")

    run_ini_text = update_run_ini_text(
        run_template_text,
        eft_ini_path.name,
        highl_template_path,
    )
    run_ini_path = run_dir / f"params_{sample_id:05d}.ini"
    run_ini_path.write_text(run_ini_text, encoding="utf-8")

    return file_root, run_ini_path, eft_ini_path


def write_run_ini_inplace(params, sample_id, run_ini_archive_dir, save_run_ini):
    """
    原始模式：覆盖写回 params_EFT.ini，按需归档一份。
    """
    if not TEMPLATE_EFT_INI.is_file():
        raise FileNotFoundError(f"找不到模板 EFT ini: {TEMPLATE_EFT_INI}")

    text = TEMPLATE_EFT_INI.read_text(encoding="utf-8")

    file_root = f"hsamp_{sample_id:05d}"
    new_text = update_ini_text(text, params, file_root)

    TEMPLATE_EFT_INI.write_text(new_text, encoding="utf-8")

    if save_run_ini and run_ini_archive_dir is not None:
        run_ini_archive_dir.mkdir(parents=True, exist_ok=True)
        archive_path = run_ini_archive_dir / f"params_EFT_{sample_id:05d}.ini"
        archive_path.write_text(new_text, encoding="utf-8")

    return file_root


# ===== 3. 调用 ./camb 运行一次 EFTCAMB =====

def run_camb(camb_exe, run_ini_path, work_dir, omp_threads, timeout_seconds=CAMB_TIMEOUT_SECONDS):
    """
    调用 ./camb run_ini_path。
    返回 (success, returncode, stdout, stderr)。
    success 仅根据 returncode 判断，进一步的物理稳定性可以在 stdout 里再查。
    若运行超过 timeout_seconds，则视为失败并返回 returncode=124。
    """
    cmd = [str(camb_exe), str(run_ini_path)]
    env = os.environ.copy()
    if omp_threads is not None:
        env.setdefault("OMP_NUM_THREADS", str(omp_threads))

    try:
        result = subprocess.run(
            cmd,
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            env=env,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        timeout_msg = f"CAMB timed out after {timeout_seconds} seconds"
        if stderr:
            stderr = f"{timeout_msg}\n{stderr}"
        else:
            stderr = timeout_msg
        return False, 124, stdout, stderr

    success = (result.returncode == 0)
    return success, result.returncode, result.stdout, result.stderr


# ===== 4. 根据工作目录找到背景输出并保存为 Horndeski_sample_i.dat =====

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

def collect_background(run_dir, sample_id, samples_dir):
    """
    从固定路径 Horndeski_solution.dat 复制为 Horndeski_sample_{sample_id}.dat
    """
    samples_dir.mkdir(parents=True, exist_ok=True)

    # 固定输出文件
    background_path = run_dir / "Horndeski_solution.dat"

    if not background_path.is_file():
        raise FileNotFoundError(f"找不到背景输出文件: {background_path}")

    target_path = samples_dir / f"Horndeski_sample_{sample_id}.dat"
    shutil.copy2(background_path, target_path)
    return target_path


def run_sample(
    sample_id,
    sample_index,
    params,
    omega_m_draw,
    run_ini_archive_dir,
    save_run_ini,
    samples_dir,
    run_work_dir_base,
    eft_template_text,
    run_template_text,
    highl_template_path,
    camb_exe,
    omp_threads,
    sample_log_path,
    sample_log_lock,
):
    log_lines = []
    log_lines.append(f"\n=== Sample {sample_index}/{N_SAMPLES} ===")
    log_lines.append(f"  Omega_m_draw: {omega_m_draw}")
    log_lines.append(f"  parameters: {params}")

    run_dir = make_run_dir(run_work_dir_base, sample_id)
    log_lines.append(f"  run dir: {run_dir}")

    file_root, run_ini_path, eft_ini_path = write_run_files(
        params,
        sample_id,
        run_dir,
        run_ini_archive_dir,
        save_run_ini,
        eft_template_text,
        run_template_text,
        highl_template_path,
    )
    log_lines.append(f"  file_root: {file_root}")

    ini_paths = (run_ini_path, eft_ini_path)
    try:
        success, returncode, out, err = run_camb(
            camb_exe,
            run_ini_path,
            run_dir,
            omp_threads,
        )
        log_lines.append(f"  [DEBUG] returncode: {returncode}")

        if not success:
            log_lines.append("  CAMB run failed, returncode != 0")
            log_lines.append("  stderr (first lines):")
            for line in err.splitlines()[:5]:
                log_lines.append(f"    {line}")
            append_sample_log(
                sample_log_path,
                sample_log_lock,
                sample_id,
                params,
                omega_m_draw,
                False,
            )
            return {"success": False, "log_lines": log_lines}

        log_lines.append("  CAMB stdout (first lines):")
        for line in out.splitlines()[:5]:
            log_lines.append(f"    {line}")

        try:
            target_path = collect_background(run_dir, sample_id, samples_dir)
            log_lines.append(f"  background saved to {target_path}")
            append_sample_log(
                sample_log_path,
                sample_log_lock,
                sample_id,
                params,
                omega_m_draw,
                True,
            )
            return {"success": True, "log_lines": log_lines}
        except FileNotFoundError as e:
            log_lines.append(f"  WARNING: {e}")
            append_sample_log(
                sample_log_path,
                sample_log_lock,
                sample_id,
                params,
                omega_m_draw,
                False,
            )
            return {"success": False, "log_lines": log_lines}
    finally:
        if not save_run_ini:
            for ini_path in ini_paths:
                try:
                    ini_path.unlink()
                except FileNotFoundError:
                    pass
                except OSError as exc:
                    log_lines.append(f"  WARNING: failed to remove ini {ini_path}: {exc}")



# ===== 5. 主循环：采样 → 写 ini → 跑 camb → 检查 → 收集 =====

def main():
    args = parse_args()
    use_parallel = resolve_bool(args.parallel, "USE_PARALLEL", False)
    save_run_ini = resolve_bool(args.save_run_ini, "SAVE_RUN_INI", False)

    if save_run_ini:
        run_ini_archive_dir = resolve_dir(
            args.run_ini_archive_dir,
            "RUN_INI_ARCHIVE_DIR",
            DEFAULT_RUN_INI_ARCHIVE_DIR,
        )
    else:
        run_ini_archive_dir = None
    samples_dir = resolve_dir(
        args.samples_dir,
        "SAMPLES_DIR",
        DEFAULT_SAMPLES_DIR,
    )

    if use_parallel:
        run_work_dir_base = resolve_dir(
            args.run_work_dir_base,
            "RUN_WORK_DIR_BASE",
            DEFAULT_RUN_WORK_DIR_BASE,
        )
        default_workers = os.cpu_count() or 1
    else:
        run_work_dir_base = None
        default_workers = 1
    max_workers = resolve_int(
        args.max_workers,
        "MAX_WORKERS",
        default_workers,
    )
    if max_workers < 1:
        raise ValueError("max_workers 必须 >= 1")
    if not use_parallel:
        max_workers = 1

    omp_threads = args.omp_threads
    if omp_threads is None and use_parallel and max_workers > 1:
        omp_threads = 1

    if not PARAM_RANGES:
        raise RuntimeError("PARAM_RANGES 为空，请先在脚本顶部填入 Horndeski 参数名及采样范围。")

    if not CAMB_EXE.is_file():
        raise FileNotFoundError(f"找不到可执行文件: {CAMB_EXE}")

    if not TEMPLATE_EFT_INI.is_file():
        raise FileNotFoundError(f"找不到模板 EFT ini: {TEMPLATE_EFT_INI}")

    if not TEMPLATE_RUN_INI.is_file():
        raise FileNotFoundError(f"找不到模板运行 ini: {TEMPLATE_RUN_INI}")

    print(f"Working dir: {CAMB_DIR}")
    print(f"Template EFT ini: {TEMPLATE_EFT_INI}")
    print(f"Template run ini: {TEMPLATE_RUN_INI}")
    if save_run_ini:
        print(f"Run ini archive dir: {run_ini_archive_dir}")
    else:
        print("Run ini archive dir: disabled")
    print(f"Samples dir: {samples_dir}")
    print(f"Parallel enabled: {use_parallel}")
    if use_parallel:
        print(f"Run work dir base: {run_work_dir_base}")
    print(f"Max workers: {max_workers}")
    print(f"OMP threads: {omp_threads}")
    print(f"Total samples: {N_SAMPLES}")
    sample_id_start = resolve_int(args.sample_id_start, "SAMPLE_ID_START", 1)
    if sample_id_start < 1:
        raise ValueError("sample_id_start 必须 >= 1")
    print(f"Sample id start: {sample_id_start}")
    readme_path, run_id = append_run_readme(samples_dir, sample_id_start)
    print(f"Run README: {readme_path}")
    sample_log_path = init_sample_log(samples_dir)
    print(f"Sample log: {sample_log_path}")

    success_count = 0
    processed_count = 0
    def update_readme(status):
        update_run_readme_progress(
            samples_dir,
            run_id,
            success_count,
            processed_count,
            status,
        )

    if use_parallel:
        highl_template_path = CAMB_DIR / "HighLExtrapTemplate_lenspotentialCls.dat"
        if not highl_template_path.is_file():
            raise FileNotFoundError(f"找不到 HighL 模板文件: {highl_template_path}")

        eft_template_text = TEMPLATE_EFT_INI.read_text(encoding="utf-8")
        run_template_text = TEMPLATE_RUN_INI.read_text(encoding="utf-8")
        sample_log_lock = threading.Lock()

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            pending = set()
            future_to_sample = {}
            try:
                for offset in range(N_SAMPLES):
                    sample_id = sample_id_start + offset
                    sample_index = offset + 1
                    params, omega_m_draw = sample_parameters()
                    future = executor.submit(
                        run_sample,
                        sample_id,
                        sample_index,
                        params,
                        omega_m_draw,
                        run_ini_archive_dir,
                        save_run_ini,
                        samples_dir,
                        run_work_dir_base,
                        eft_template_text,
                        run_template_text,
                        highl_template_path,
                        CAMB_EXE,
                        omp_threads,
                        sample_log_path,
                        sample_log_lock,
                    )
                    pending.add(future)
                    future_to_sample[future] = (sample_id, sample_index)

                    if len(pending) >= max_workers * 2:
                        done, pending = wait(pending, return_when=FIRST_COMPLETED)
                        for fut in done:
                            sample_id, sample_index = future_to_sample.pop(fut)
                            try:
                                result = fut.result()
                            except Exception as e:
                                print(f"\n=== Sample {sample_index}/{N_SAMPLES} ===")
                                print(f"  ERROR: {e}")
                                processed_count += 1
                                update_readme("running")
                                continue
                            print("\n".join(result["log_lines"]))
                            if result["success"]:
                                success_count += 1
                            processed_count += 1
                            update_readme("running")

                while pending:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)
                    for fut in done:
                        sample_id, sample_index = future_to_sample.pop(fut)
                        try:
                            result = fut.result()
                        except Exception as e:
                            print(f"\n=== Sample {sample_index}/{N_SAMPLES} ===")
                            print(f"  ERROR: {e}")
                            processed_count += 1
                            update_readme("running")
                            continue
                        print("\n".join(result["log_lines"]))
                        if result["success"]:
                            success_count += 1
                        processed_count += 1
                        update_readme("running")
            except KeyboardInterrupt:
                print("\nInterrupted by user.")
                update_readme("interrupted")
                executor.shutdown(wait=False, cancel_futures=True)
                return
    else:
        try:
            for offset in range(N_SAMPLES):
                sample_id = sample_id_start + offset
                sample_index = offset + 1
                print(f"\n=== Sample {sample_index}/{N_SAMPLES} ===")

                params, omega_m_draw = sample_parameters()
                print("  Omega_m_draw:", omega_m_draw)
                print("  parameters:", params)

                file_root = write_run_ini_inplace(
                    params,
                    sample_id,
                    run_ini_archive_dir,
                    save_run_ini,
                )
                print("  file_root:", file_root)

                success, returncode, out, err = run_camb(
                    CAMB_EXE,
                    TEMPLATE_RUN_INI,
                    CAMB_DIR,
                    omp_threads,
                )
                print("  [DEBUG] returncode:", returncode)
                for line in out.splitlines()[:10]:
                    print("    [STDOUT]", line)

                if not success:
                    print("  CAMB run failed, returncode != 0")
                    print("  stderr (first lines):")
                    for line in err.splitlines()[:5]:
                        print("    ", line)
                    append_sample_log(
                        sample_log_path,
                        None,
                        sample_id,
                        params,
                        omega_m_draw,
                        False,
                    )
                    processed_count += 1
                    update_readme("running")
                    continue

                print("  CAMB stdout (first lines):")
                for line in out.splitlines()[:5]:
                    print("    ", line)

                try:
                    background_path = CAMB_DIR / "Horndeski_solution.dat"
                    stat = background_path.stat()
                    print("  [DEBUG] Horndeski_solution.dat size:", stat.st_size)
                    target_path = collect_background(CAMB_DIR, sample_id, samples_dir)
                    print(f"  background saved to {target_path}")
                    success_count += 1
                    append_sample_log(
                        sample_log_path,
                        None,
                        sample_id,
                        params,
                        omega_m_draw,
                        True,
                    )
                except FileNotFoundError as e:
                    print("  WARNING:", e)
                    append_sample_log(
                        sample_log_path,
                        None,
                        sample_id,
                        params,
                        omega_m_draw,
                        False,
                    )
                processed_count += 1
                update_readme("running")
        except KeyboardInterrupt:
            print("\nInterrupted by user.")
            update_readme("interrupted")
            return

    update_readme("complete")
    print(f"\nDone. Accepted {success_count} / {N_SAMPLES} samples.")


if __name__ == "__main__":
    main()
