# -*- coding: utf-8 -*-
# Author; alin m elena, alin@elena.re
# Contribs;
# Date: 01-02-2024
# ©alin m elena, GPL v3 https://www.gnu.org/licenses/gpl-3.0.en.html


from pathlib import Path

from urllib.parse import urlparse
from urllib import request
import argparse

cli=argparse.ArgumentParser()
cli.add_argument(
  "--url",
  type=str,
  default='http://tinyurl.com/46jrkm3v',
  help = ' download a model from the url. default %(default)s'
  )

cli.add_argument(
  "--path",
  type=str,
  default='~/.cache/mace',
  help = ' save the model in the folder. default %(default)s'
  )

args = cli.parse_args()

url = args.url
path = args.path

p = Path(path).expanduser()
model_name = urlparse(url).path.split("/")[-1]
p.mkdir(parents=True, exist_ok=True)

save_path_model = p.joinpath(model_name)

if save_path_model.exists():
  print(f"{model_name} already cached in {save_path_model}")
else:
  h, msg = request.urlretrieve(url, save_path_model)
  print(msg)
  if "Content-Type" in msg:
    raise RuntimeError(f"wrong content download from {url} check the model link")
  else:
    print(f"download model {model_name} from {url} and save it in {save_path_model}")
print(f"all cached models in {p}")
for f in p.glob("*"):
    print(f)
