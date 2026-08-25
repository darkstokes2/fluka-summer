# FLUKA Batch Processor (USRBIN)

批处理 FLUKA USRBIN 模拟输出：用 `usbsuw` 合并 `_fort.NN` 文件，用 `usbrea` 转成 `.lis`。

## 功能

- 自动发现或手动指定模拟文件夹
- pexpect 驱动 `usbsuw` 合并多个 `_fort.NN` → 单个 `.bnn`
- pexpect 驱动 `usbrea` 转换 `.bnn` → `.lis`
- 全部命令行参数可配置，无硬编码路径

## 依赖

```
pexpect
```

还需要 FLUKA 工具链：`usbsuw`、`usbrea`（需在 PATH 中或通过参数指定路径）。

## 用法

```bash
# 自动发现 base 目录下含 _fort.42 的子文件夹
python fluka_batch_processor.py --base /path/to/simulation/projects

# 指定文件夹列表和 fort 单元号
python fluka_batch_processor.py --base /path/to/projects \
    --folders run_A run_B --fort-unit 42 --tag Proton

# 指定 FLUKA 工具路径和超时
python fluka_batch_processor.py --base /path/to/projects \
    --usbsuw-bin /opt/fluka/bin/usbsuw \
    --usbrea-bin /opt/fluka/bin/usbrea \
    --usbsuw-timeout 600 --usbrea-timeout 1200
```

### 参数说明

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `--base` | （必填） | 包含模拟子文件夹的基础目录 |
| `--post` | `<base>/postprocess` | `.lis` 输出目录 |
| `--folders` | 自动发现 | 要处理的子文件夹名列表 |
| `--fort-unit` | `42` | Fortran 单元号（42 = 能量沉积） |
| `--tag` | `Proton` | 输出文件命名标签 |
| `--usbsuw-bin` | `usbsuw` | usbsuw 可执行文件路径 |
| `--usbrea-bin` | `usbrea` | usbrea 可执行文件路径 |
| `--usbsuw-timeout` | `300` | usbsuw 超时（秒） |
| `--usbrea-timeout` | `600` | usbrea 超时（秒） |

### 输出

- `.bnn` 文件留在各自源文件夹
- `.lis` 文件输出到 `--post` 目录

## 处理流程

```
对每个文件夹：
  1. 查找 *_fort.NN 文件（无则跳过）
  2. cd 到目标文件夹
  3. usbsuw 合并 → .bnn
  4. usbrea 转换 → .lis（输出到 post 目录）
  5. 记录成功/失败状态
```

生成的 `.lis` 文件可直接供 `usrbin-reader` 使用。
