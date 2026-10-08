# By serving endpoint

Printed by `analyze.py` (WS-C1) from the default configuration (`scores/main`), the table of RESULTS.md's "By serving endpoint" (D-133, D-149) with the templates each endpoint served. Open-weight models were served by several endpoints under the fp8-or-better rule; the harness dispatches items in template order and OpenRouter falls back under load, so a raw per-endpoint mean is confounded with the templates each endpoint happened to serve. The matched difference compares an endpoint with the other endpoints of the same model on the same templates. Two flags: the endpoint served fewer than 20 templates; its matched difference rests on fewer than 20 templates, which also happens to the endpoint that served nearly every row when the others served few templates. Models served by one endpoint are not listed.

| model | endpoint | rows | raw score | unusable | matched difference | templates matched | templates served | served fewer than 20 | matched on fewer than 20 |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| `gpt-oss-20b` | Darkbloom | 2240 | 0.808 | 0.026 | -0.228 | 9 | 150 |  | yes |
| `gpt-oss-20b` | DekaLLM | 9 | 0.889 | 0.000 | +0.239 | 8 | 8 | yes | yes |
| `gpt-oss-20b` | DeepInfra | 1 | 1.000 | 0.000 | +0.143 | 1 | 1 | yes | yes |
| `gemma-4-26b-a4b` | DekaLLM | 1245 | 0.883 | 0.000 | +0.002 | 148 | 148 |  |  |
| `gemma-4-26b-a4b` | NextBit | 975 | 0.875 | 0.000 | -0.002 | 148 | 148 |  |  |
| `gemma-4-26b-a4b` | Io Net | 30 | 0.333 | 0.000 |  | 0 | 2 | yes | yes |
| `deepseek-v4.1-flash` | CoreWeave | 782 | 0.992 | 0.001 | -0.000 | 141 | 141 |  |  |
| `deepseek-v4.1-flash` | Parasail | 466 | 0.992 | 0.000 | -0.002 | 138 | 138 |  |  |
| `deepseek-v4.1-flash` | Morph | 456 | 0.991 | 0.000 | +0.002 | 133 | 133 |  |  |
| `deepseek-v4.1-flash` | DeepInfra | 411 | 0.984 | 0.005 | +0.000 | 78 | 78 |  |  |
| `deepseek-v4.1-flash` | Makora | 87 | 1.000 | 0.000 | +0.003 | 43 | 43 |  |  |
| `deepseek-v4.1-flash` | InferenceNet | 30 | 0.800 | 0.167 |  | 0 | 2 | yes | yes |
| `deepseek-v4.1-flash` | Novita | 16 | 1.000 | 0.000 | +0.016 | 7 | 7 | yes | yes |
| `deepseek-v4.1-flash` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 | 2 | yes | yes |
| `qwen3-235b-a22b-2507` | GMICloud | 2192 | 0.899 | 0.000 | -0.010 | 35 | 150 |  |  |
| `qwen3-235b-a22b-2507` | Parasail | 54 | 0.907 | 0.000 | +0.006 | 35 | 35 |  |  |
| `qwen3-235b-a22b-2507` | Nebius | 4 | 1.000 | 0.000 | +0.182 | 1 | 1 | yes | yes |
| `glm-5.3-flash` | GMICloud | 981 | 0.976 | 0.015 | -0.003 | 145 | 145 |  |  |
| `glm-5.3-flash` | Novita | 419 | 0.962 | 0.031 | +0.003 | 120 | 121 |  |  |
| `glm-5.3-flash` | Parasail | 220 | 0.964 | 0.036 | +0.002 | 75 | 75 |  |  |
| `glm-5.3-flash` | Sail Research | 201 | 0.965 | 0.020 | -0.022 | 114 | 114 |  |  |
| `glm-5.3-flash` | Phala | 166 | 0.991 | 0.006 | +0.004 | 74 | 74 |  |  |
| `glm-5.3-flash` | StreamLake | 80 | 0.950 | 0.037 | -0.025 | 35 | 35 |  |  |
| `glm-5.3-flash` | Z.AI | 79 | 0.994 | 0.000 | +0.004 | 48 | 48 |  |  |
| `glm-5.3-flash` | NextBit | 43 | 0.965 | 0.023 | +0.031 | 13 | 13 | yes | yes |
| `glm-5.3-flash` | AtlasCloud | 42 | 1.000 | 0.000 | +0.014 | 15 | 15 | yes | yes |
| `glm-5.3-flash` | SiliconFlow | 18 | 1.000 | 0.000 | +0.000 | 5 | 5 | yes | yes |
| `glm-5.3-flash` | Morph | 1 | 1.000 | 0.000 | +0.000 | 1 | 1 | yes | yes |
| `glm-5.3` | Baidu | 918 | 0.956 | 0.036 | -0.006 | 86 | 95 |  |  |
| `glm-5.3` | Sail Research | 739 | 0.964 | 0.030 | -0.004 | 77 | 80 |  |  |
| `glm-5.3` | Morph | 257 | 0.940 | 0.054 | +0.021 | 123 | 123 |  |  |
| `glm-5.3` | Novita | 167 | 0.943 | 0.054 | -0.004 | 41 | 41 |  |  |
| `glm-5.3` | AtlasCloud | 137 | 0.927 | 0.073 | -0.010 | 54 | 54 |  |  |
| `glm-5.3` | AkashML | 15 | 0.267 | 0.733 | -0.020 | 5 | 5 | yes | yes |
| `glm-5.3` | BaseTen | 15 | 1.000 | 0.000 | +0.000 | 5 | 5 | yes | yes |
| `glm-5.3` | GMICloud | 2 | 1.000 | 0.000 | +0.000 | 2 | 2 | yes | yes |
| `kimi-k3` | Morph | 2231 | 0.980 | 0.007 | -0.017 | 7 | 150 |  | yes |
| `kimi-k3` | BaseTen | 19 | 1.000 | 0.000 | +0.017 | 7 | 7 | yes | yes |
| `gpt-5.4-mini` | Azure | 2249 | 0.852 | 0.000 | -0.036 | 1 | 150 |  | yes |
| `gpt-5.4-mini` | OpenAI | 1 | 1.000 | 0.000 | +0.036 | 1 | 1 | yes | yes |

Over the 21 endpoints whose matched difference rests on at least 20 templates, it lies between -0.025 and +0.021 (largest absolute value 0.025); of the other 19, 16 served fewer than 20 templates.
