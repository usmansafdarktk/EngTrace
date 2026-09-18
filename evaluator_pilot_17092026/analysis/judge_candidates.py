"""Judge candidates for E1/E5: every OpenRouter model outside the evaluated families.

The evaluated suite is the paper's Table 13 (27 models), which spans seven
families once backbones are counted: OpenAI, Anthropic, Google (Gemini AND Gemma),
DeepSeek, Meta (Llama, and MetaMath via Llemma), Alibaba/Qwen, and Mistral
(Mathstral, and WizardMath via Mistral-7B). A judge from any of them reintroduces
the judge/judged overlap E1 exists to remove.

For each remaining family this prints what can be read off the catalogue: current
ids, price, context, whether weights are published (OpenRouter's hugging_face_id),
and whether structured output is supported - the Tribunal requires a JSON object.
Capability and lineage are not in the catalogue and are assessed separately.
"""
import os as _os
_ANALYSIS = _os.path.dirname(_os.path.abspath(__file__))
_REPO = _os.path.dirname(_os.path.dirname(_ANALYSIS))

import datetime as dt
import json
import os
import urllib.request

from dotenv import load_dotenv

load_dotenv(_os.path.join(_REPO, '.env'))
EXCLUDED = ('openai/', 'anthropic/', 'google/', 'deepseek/', 'meta-llama/', 'qwen/',
            'mistralai/', 'alibaba/', 'tngtech/')          # tng: deepseek merges
H = {'Authorization': 'Bearer ' + os.environ['OPENROUTER_API_KEY'], 'User-Agent': 'engtrace'}
cat = json.load(urllib.request.urlopen(
    urllib.request.Request('https://openrouter.ai/api/v1/models', headers=H), timeout=60))['data']

by_org = {}
for m in cat:
    mid = m['id']
    if ':' in mid or mid.startswith('~') or mid.startswith(EXCLUDED) or '/' not in mid:
        continue
    by_org.setdefault(mid.split('/')[0], []).append(m)

print('Families outside the evaluated suite, newest 4 text models each:\n')
for org in sorted(by_org):
    ms = sorted(by_org[org], key=lambda m: -(m.get('created') or 0))
    ms = [m for m in ms if 'text' in (m.get('architecture', {}).get('output_modalities') or ['text'])][:4]
    if not ms:
        continue
    print(org)
    for m in ms:
        p = m['pricing']
        sp = set(m.get('supported_parameters') or [])
        print('  %-44s %s  in $%5.2f out $%6.2f  ctx %7s  %s  %s' % (
            m['id'][:44],
            dt.datetime.fromtimestamp(m.get('created') or 0, dt.timezone.utc).strftime('%Y-%m'),
            float(p.get('prompt') or 0) * 1e6, float(p.get('completion') or 0) * 1e6,
            m.get('context_length'),
            'open' if m.get('hugging_face_id') else 'closed',
            'json' if ({'response_format', 'structured_outputs'} & sp) else '-'))
