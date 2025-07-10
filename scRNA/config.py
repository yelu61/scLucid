import yaml
import os
from typing import Dict, Any

# 默认配置
DEFAULT_CONFIG = {
    "qc": {
        "min_genes": 200,
        "max_genes": 6000,
        "pc_mt": 20,
        "pc_hb": 20
    },
    "norm": {
        "target_sum": 1e4,
        "exclude_highly_expressed": False
    },
    "hvg": {
        "method": "scanpy",
        "flavor": "seurat",
        "n_top_genes": 2000
    },
    # 其他默认配置...
}

def load_config(config_path: str = None) -> Dict[str, Any]:
    """加载配置文件，如不存在则使用默认配置"""
    config = DEFAULT_CONFIG.copy()
    
    if config_path and os.path.exists(config_path):
        with open(config_path, 'r') as f:
            user_config = yaml.safe_load(f)
            # 递归更新配置
            _update_config(config, user_config)
    
    return config

def _update_config(default_config, user_config):
    """递归更新配置字典"""
    for key, value in user_config.items():
        if isinstance(value, dict) and key in default_config:
            _update_config(default_config[key], value)
        else:
            default_config[key] = value