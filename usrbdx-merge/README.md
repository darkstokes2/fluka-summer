# USRBDX Batch Processor

批处理 FLUKA USRBDX 模拟输出：用 `usxsuw` 合并 `_fort.NN` 文件，自动生成 `.bnn` + `_sum.lis` + `_tab.lis`。

## 功能

- 自动发现或手动指定模拟文件夹
- 支持多探测器并行处理（`--fort-units 50 51 52`）
- pexpect 驱动 `usxsuw`，一步生成三个输出文件（无需 usbrea）
- 全部命令行参数可配置，无硬编码路径

## 依赖

```
pexpect
```

还需要 FLUKA 工具 `usxsuw`（需在 PATH 中或通过参数指定路径）。

## 用法

```bash
# 自动发现 base 目录下含 _fort.50 的子文件夹，处理 Proton 探测器
python usrbdx_batch_processor.py --base /path/to/simulation/projects

# 指定文件夹列表，处理多个探测器
python usrbdx_batch_processor.py --base /path/to/projects \
    --folders run_A run_B \
    --fort-units 50 51 52 --tags Proton Neutron Electron

# 指定工具路径和超时
python usrbdx_batch_processor.py --base /path/to/projects \
    --usxsuw-bin /opt/fluka/bin/usxsuw --usxsuw-timeout 600
```

### 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--base` | （必填） | 包含模拟子文件夹的基础目录 |
| `--post` | `<base>/postprocess` | `.lis` 输出目录 |
| `--folders` | 自动发现 | 要处理的子文件夹名列表 |
| `--fort-units` | `50` | USRBDX Fortran 单元号列表 |
| `--tags` | `Proton` | 对应标签列表（与 --fort-units 等长） |
| `--usxsuw-bin` | `usxsuw` | usxsuw 可执行文件路径 |
| `--usxsuw-timeout` | `300` | usxsuw 超时（秒） |

### 输出

- `.bnn` 文件留在各自源文件夹
- `_sum.lis` 和 `_tab.lis` 移到 `--post` 目录

### 与 USRBIN 批处理的区别

| 方面 | USRBIN 批处理 | USRBDX 批处理 |
|------|-------------|-------------|
| FLUKA 工具 | `usbsuw` + `usbrea`（两步） | `usxsuw`（一步） |
| 输入 | `_fort.42` | `_fort.50/51/52` |
| 输出 | `.bnn` + `.lis` | `.bnn` + `_sum.lis` + `_tab.lis` |
| 每文件夹探测器数 | 单单元 | 多单元（并行列表） |

生成的 `_tab.lis` 文件可直接供 `usrbdx-reader` 使用。
