# output_diff -- 输出级逐字节回归

> 只收录**需要输出级逐字节比对**的用例；单元/CLI 断言（细粒度）见 `../unit`、`../cli`。
> 设计借鉴 `us-align_modify/.../cLanguage2Cplus` 的 baseline/current/diffs 框架，适配 LS-align
> "部分改动必须逐字节一致、部分改动（R 批）是故意改数值"的双情形。

## 定位

| 目录 | 职责 |
|---|---|
| `../unit/` | 白盒单测（46 项），函数级断言 |
| `../cli/` | 黑盒 CLI（42 项），行为/退出码/表格结构 |
| `../golden/` | 2 个刚性快照（`ace_rigid`、`self_copy`），"必须一致"的最小集 |
| **`output_diff/`** | **扩展的完整输出快照回归**：多用例、分 `exact` / `expected` 两模式、自动剥离时序行、生成 diff |

## 结构

```
output_diff/
├── cases.txt          # 用例清单：name | workdir | args | mode
├── output_diff.py     # runner（Python；剥离时序行 + 归一化 + 分类比对）
├── baseline/          # 黄金标准输出（提交进 git，不可随意改）
├── current/           # 每次运行输出（可删）
├── diffs/             # 不匹配时的 unified diff（可删）
├── expected/          # R 批次"期望差异清单"（JSON）
└── README.md
```

## 用例模式（`cases.txt` 末列）

| mode | 含义 | 判定 |
|---|---|---|
| `exact` | 业务数据必须与 baseline 逐字节一致（A 批 / F 系列 / 等价重构） | 任何业务差异 = **FAIL** |
| `expected:KEY` | 输出**允许变化**，但须符合 `expected/KEY.json`（R 批故意改数值） | 无变化 = FAIL（说明没生效）；有变化 = WARN（留 diff 人工核对） |

## 用法

```bash
cd tests

# 1) 建立/重建基线（在"已知正确"的提交上运行一次，提交 baseline/）
python output_diff/output_diff.py --baseline

# 2) 改动后比对
python output_diff/output_diff.py

# 3) 只跑一个用例
python output_diff/output_diff.py --case ace_rigid
```

退出码：全 PASS → 0；有 FAIL → 1。

## 归一化（比对前自动做）

- **剥离时序行**：`Using 0.003-s to load mol2 files.` / `Using 0.001-s to do LS-align.` / `Totally Using 0.004-s.`（数值每次不同）
- **换行归一化**：`\r\n` / 单独 `\r` / 连续换行统一为 `\n`
- **写入固定 LF**：`current/`、`baseline/`、`diffs/` 均以 `newline="\n"` 写出，Windows 与 Linux 跑出的 baseline 字节一致（与 `tests/**` 的 LF 约定相同）
- 路径前缀差异（workdir 不同）不体现在表体里，无需额外处理

## 与 R 批（会改数值）的配合

1. **修 R 批前**：相关用例保持 `exact`，baseline 是修复前的输出。
2. **修 R 批后**：该用例的输出**应当变化**——因此：
   - 若要保留"变化符合预期"的自动判定 → 改成 `expected:<key>` + 写 `expected/<key>.json`；
   - 若某项修复**不改 stdout**（如 R8 只改搜索矩阵、报告分不变）→ 保持 `exact`，正确性由 `unit/` 的 `r8_*` 用例保证，不硬塞进 output_diff。
3. **重生成 baseline**：确认新输出符合设计后，`--baseline` 重建并提交，作为下一轮基线。

> 重要：baseline 必须由**目标提交**真实运行产生，不能用"改后程序"去覆盖旧基线后声称一致。

## 与 golden 的关系

`golden/`（`make check-golden`）是 output_diff 的**子集**（2 个刚性快照），保留它是因为它已接入既有 CI（`make test`）。
output_diff 提供更广的用例、模式区分和 diff 归档；两者不冲突，`golden` 不变即 output_diff 对应两项也应 PASS。

## 环境注意

- 二进制在 `../build/LSalign_under_test`（`make` 产物）。
- LS-align **始终退出 1**（成功也返回 1），本 runner 不依赖退出码判定成败，只凭输出比对。
- 柔性用例（`-rf 1`）含 `srand(clock())` 随机抽样——`octane_flex`/`ace_flex` 之所以能纳入 `exact`，是因为当前数据下抽样路径确定（构象数 < 5000，不触发随机子集）。若将来用例触发随机抽样，应移出 `exact` 或固定种子。
