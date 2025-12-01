"""
Global Storage Manager for scLucid Package

Provides unified management of analysis results stored in adata.uns['sclucid'],
including status checking, cleanup, and optimization functions.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union
import pandas as pd
from pathlib import Path
import pickle
from datetime import datetime
from anndata import AnnData

log = logging.getLogger(__name__)

class StorageManager:
    """Global storage manager for scLucid analysis results"""
    
    # 定义模块和存储路径的映射
    MODULE_STORAGE_MAP = {
        'de': 'sclucid.analysis.de',  # 差异分析
        'enrichment': 'sclucid.analysis.de.enrichment',  # 富集分析
        'clustering': 'sclucid.analysis.clustering',  # 聚类分析
        'annotation': 'sclucid.analysis.annotation',  # 细胞注释
        'scoring': 'sclucid.analysis.scoring',  # 基因集评分
        'proportion': 'sclucid.analysis.proportion',  # 细胞比例分析
        'visualization': 'sclucid.visualization',  # 可视化配置
        'characterization': 'cluster_characterization',  # 簇特征化
    }
    
    # 定义结果类型和对应的键模式
    RESULT_TYPE_PATTERNS = {
        'raw_results': ['_raw$', '_original$'],  # 原始结果
        'processed_results': ['_df$', '_filtered$', '_processed$'],  # 处理后结果
        'config': ['_config$', '_params$'],  # 配置参数
        'integrated': ['characterization', 'summary'],  # 整合结果
        'temp': ['^temp_', '^_temp_'],  # 临时结果
    }
    
    @classmethod
    def get_storage_status(
        cls, 
        adata: AnnData, 
        module: Optional[str] = None,
        include_size: bool = True,
        include_details: bool = False
    ) -> Dict:
        """
        获取存储状态概览
        
        Parameters:
            adata: AnnData对象
            module: 指定模块，None表示所有模块
            include_size: 是否计算存储大小
            include_details: 是否包含详细信息
            
        Returns:
            存储状态字典
        """
        status = {
            'total_modules': 0,
            'modules': {},
            'summary': {
                'total_keys': 0,
                'total_size_mb': 0,
                'duplicate_groups': [],
                'temp_keys': []
            }
        }
        
        if 'sclucid' not in adata.uns:
            return status
        
        # 遍历所有存储路径
        for mod_name, storage_path in cls.MODULE_STORAGE_MAP.items():
            if module and mod_name != module:
                continue
                
            mod_status = cls._get_module_status(
                adata, mod_name, storage_path, include_size, include_details
            )
            
            if mod_status['keys']:
                status['modules'][mod_name] = mod_status
                status['total_modules'] += 1
                status['summary']['total_keys'] += len(mod_status['keys'])
                status['summary']['total_size_mb'] += mod_status.get('size_mb', 0)
                
                # 收集重复组和临时键
                status['summary']['duplicate_groups'].extend(mod_status.get('duplicate_groups', []))
                status['summary']['temp_keys'].extend(mod_status.get('temp_keys', []))
        
        # 去重
        status['summary']['duplicate_groups'] = list(set(status['summary']['duplicate_groups']))
        status['summary']['temp_keys'] = list(set(status['summary']['temp_keys']))
        
        return status
    
    @classmethod
    def _get_module_status(
        cls, 
        adata: AnnData, 
        module_name: str, 
        storage_path: str,
        include_size: bool,
        include_details: bool
    ) -> Dict:
        """获取单个模块的存储状态"""
        module_status = {
            'name': module_name,
            'storage_path': storage_path,
            'keys': [],
            'size_mb': 0,
            'result_types': {},
            'duplicate_groups': [],
            'temp_keys': []
        }
        
        # 获取存储的数据
        storage_data = cls._get_nested_data(adata.uns, storage_path)
        if not storage_data:
            return module_status
        
        # 分析每个键
        for key, value in storage_data.items():
            module_status['keys'].append(key)
            
            # 计算大小
            if include_size:
                size_mb = cls._calculate_size(value)
                module_status['size_mb'] += size_mb
            
            # 分类结果类型
            result_type = cls._classify_result_type(key)
            if result_type not in module_status['result_types']:
                module_status['result_types'][result_type] = []
            module_status['result_types'][result_type].append(key)
            
            # 检查重复组
            if cls._is_duplicate_group(key):
                module_status['duplicate_groups'].append(key)
            
            # 检查临时键
            if cls._is_temp_key(key):
                module_status['temp_keys'].append(key)
            
            # 详细信息
            if include_details:
                module_status['details'] = module_status.get('details', {})
                module_status['details'][key] = cls._get_value_details(value)
        
        return module_status
    
    @classmethod
    def cleanup_storage(
        cls,
        adata: AnnData,
        module: Optional[str] = None,
        result_types: Optional[List[str]] = None,
        keys_to_keep: Optional[List[str]] = None,
        dry_run: bool = False,
        force: bool = False
    ) -> Dict:
        """
        清理存储空间
        
        Parameters:
            adata: AnnData对象
            module: 指定模块，None表示所有模块
            result_types: 要清理的结果类型，如['temp', 'raw_results']
            keys_to_keep: 要保留的键，优先级高于result_types
            dry_run: 是否只显示将要删除的内容，不实际删除
            force: 是否强制删除（跳过确认）
            
        Returns:
            清理结果报告
        """
        if not dry_run and not force:
            # 交互式确认
            status = cls.get_storage_status(adata, module, include_size=True)
            print("Current storage status:")
            cls._print_storage_status(status)
            
            response = input("\nProceed with cleanup? (y/n): ")
            if response.lower() != 'y':
                return {'status': 'cancelled', 'message': 'Cleanup cancelled by user'}
        
        cleanup_report = {
            'deleted_keys': [],
            'freed_space_mb': 0,
            'errors': []
        }
        
        # 获取要删除的键
        keys_to_delete = cls._get_keys_for_cleanup(
            adata, module, result_types, keys_to_keep
        )
        
        # 执行删除
        for key_path in keys_to_delete:
            try:
                # 计算删除前的大小
                value = cls._get_nested_data(adata.uns, key_path)
                size_mb = cls._calculate_size(value)
                
                if not dry_run:
                    # 实际删除
                    cls._delete_nested_key(adata.uns, key_path)
                
                cleanup_report['deleted_keys'].append(key_path)
                cleanup_report['freed_space_mb'] += size_mb
                
            except Exception as e:
                cleanup_report['errors'].append(f"Error deleting {key_path}: {str(e)}")
        
        # 清理空字典
        if not dry_run:
            cls._cleanup_empty_dicts(adata.uns)
        
        cleanup_report['status'] = 'completed' if not dry_run else 'dry_run'
        return cleanup_report
    
    @classmethod
    def optimize_storage(
        cls,
        adata: AnnData,
        strategy: str = 'integrated_only',
        keep_configs: bool = True,
        dry_run: bool = False
    ) -> Dict:
        """
        优化存储策略
        
        Parameters:
            adata: AnnData对象
            strategy: 优化策略
                - 'integrated_only': 只保留整合结果
                - 'processed_only': 只保留处理后结果
                - 'configs_only': 只保留配置
                - 'minimal': 最小化存储（只保留关键结果）
            keep_configs: 是否保留配置信息
            dry_run: 是否只显示将要删除的内容
            
        Returns:
            优化结果报告
        """
        strategy_configs = {
            'integrated_only': {
                'keep_types': ['integrated'],
                'keep_patterns': ['characterization', 'summary']
            },
            'processed_only': {
                'keep_types': ['processed_results', 'integrated'],
                'keep_patterns': ['_df$', '_filtered$', '_processed$', 'characterization']
            },
            'configs_only': {
                'keep_types': ['config'],
                'keep_patterns': ['_config$', '_params$']
            },
            'minimal': {
                'keep_types': ['integrated'],
                'keep_patterns': ['characterization'],
                'keep_specific': ['rank_genes_groups']  # 保留关键原始结果
            }
        }
        
        if strategy not in strategy_configs:
            raise ValueError(f"Unknown strategy: {strategy}")
        
        config = strategy_configs[strategy]
        
        # 获取当前状态
        status = cls.get_storage_status(adata, include_size=True)
        
        # 确定要保留的键
        keys_to_keep = []
        
        # 按类型保留
        for result_type, keys in cls._get_all_keys_by_type(adata).items():
            if result_type in config['keep_types']:
                keys_to_keep.extend(keys)
        
        # 按模式保留
        import re
        for pattern in config['keep_patterns']:
            for key in cls._get_all_keys(adata):
                if re.search(pattern, key):
                    keys_to_keep.append(key)
        
        # 保留特定键
        if 'keep_specific' in config:
            keys_to_keep.extend(config['keep_specific'])
        
        # 去重
        keys_to_keep = list(set(keys_to_keep))
        
        # 执行清理
        return cls.cleanup_storage(
            adata,
            keys_to_keep=keys_to_keep,
            dry_run=dry_run,
            force=True
        )
    
    # === 辅助方法 ===
    
    @classmethod
    def _get_nested_data(cls, data: Dict, path: str):
        """获取嵌套字典中的数据"""
        keys = path.split('.')
        current = data
        for key in keys:
            if key in current:
                current = current[key]
            else:
                return None
        return current
    
    @classmethod
    def _delete_nested_key(cls, data: Dict, path: str):
        """删除嵌套字典中的键"""
        keys = path.split('.')
        current = data
        for key in keys[:-1]:
            current = current[key]
        del current[keys[-1]]
    
    @classmethod
    def _calculate_size(cls, value) -> float:
        """计算值的内存大小（MB）"""
        if isinstance(value, pd.DataFrame):
            return value.memory_usage(deep=True).sum() / 1024**2
        elif isinstance(value, dict):
            return sum(cls._calculate_size(v) for v in value.values()) / 1024**2
        else:
            return 0  # 简化处理
    
    @classmethod
    def _classify_result_type(cls, key: str) -> str:
        """分类结果类型"""
        import re
        for result_type, patterns in cls.RESULT_TYPE_PATTERNS.items():
            for pattern in patterns:
                if re.search(pattern, key):
                    return result_type
        return 'other'
    
    @classmethod
    def _is_duplicate_group(cls, key: str) -> bool:
        """检查是否为重复组"""
        # 简化实现，实际可以根据业务逻辑定制
        duplicate_indicators = ['_filtered', '_processed', '_df']
        return any(indicator in key for indicator in duplicate_indicators)
    
    @classmethod
    def _is_temp_key(cls, key: str) -> bool:
        """检查是否为临时键"""
        import re
        return bool(re.search(r'^temp_|^_temp_', key))
    
    @classmethod
    def _get_value_details(cls, value) -> Dict:
        """获取值的详细信息"""
        if isinstance(value, pd.DataFrame):
            return {
                'type': 'DataFrame',
                'shape': value.shape,
                'columns': list(value.columns)[:5],  # 只显示前5列
                'memory_mb': value.memory_usage(deep=True).sum() / 1024**2
            }
        elif isinstance(value, dict):
            return {
                'type': 'dict',
                'keys': list(value.keys())[:5],  # 只显示前5个键
                'size': len(value)
            }
        else:
            return {
                'type': type(value).__name__,
                'value': str(value)[:50]  # 只显示前50个字符
            }
    
    @classmethod
    def _get_keys_for_cleanup(
        cls,
        adata: AnnData,
        module: Optional[str],
        result_types: Optional[List[str]],
        keys_to_keep: Optional[List[str]]
    ) -> List[str]:
        """获取要清理的键列表"""
        all_keys = []
        
        # 获取所有键
        for mod_name, storage_path in cls.MODULE_STORAGE_MAP.items():
            if module and mod_name != module:
                continue
                
            storage_data = cls._get_nested_data(adata.uns, storage_path)
            if storage_data:
                for key in storage_data.keys():
                    full_path = f"{storage_path}.{key}"
                    all_keys.append(full_path)
        
        # 确定要删除的键
        if keys_to_keep:
            keys_to_delete = [k for k in all_keys if not any(keep in k for keep in keys_to_keep)]
        elif result_types:
            keys_to_delete = []
            for key_path in all_keys:
                key = key_path.split('.')[-1]
                result_type = cls._classify_result_type(key)
                if result_type in result_types:
                    keys_to_delete.append(key_path)
        else:
            keys_to_delete = all_keys
        
        return keys_to_delete
    
    @classmethod
    def _cleanup_empty_dicts(cls, data: Dict):
        """清理空字典"""
        keys_to_delete = []
        for key, value in data.items():
            if isinstance(value, dict):
                cls._cleanup_empty_dicts(value)
                if not value:
                    keys_to_delete.append(key)
        
        for key in keys_to_delete:
            del data[key]
    
    @classmethod
    def _print_storage_status(cls, status: Dict):
        """打印存储状态"""
        print(f"\n=== Storage Status ===")
        print(f"Total modules: {status['total_modules']}")
        print(f"Total keys: {status['summary']['total_keys']}")
        print(f"Total size: {status['summary']['total_size_mb']:.2f} MB")
        
        if status['summary']['duplicate_groups']:
            print(f"Duplicate groups: {len(status['summary']['duplicate_groups'])}")
        
        if status['summary']['temp_keys']:
            print(f"Temporary keys: {len(status['summary']['temp_keys'])}")
        
        for mod_name, mod_status in status['modules'].items():
            print(f"\n--- Module: {mod_name} ---")
            print(f"Keys: {len(mod_status['keys'])}")
            print(f"Size: {mod_status.get('size_mb', 0):.2f} MB")
            print(f"Result types: {list(mod_status['result_types'].keys())}")
    
    @classmethod
    def _get_all_keys_by_type(cls, adata: AnnData) -> Dict[str, List[str]]:
        """按类型获取所有键"""
        keys_by_type = {}
        
        for mod_name, storage_path in cls.MODULE_STORAGE_MAP.items():
            storage_data = cls._get_nested_data(adata.uns, storage_path)
            if storage_data:
                for key in storage_data.keys():
                    result_type = cls._classify_result_type(key)
                    if result_type not in keys_by_type:
                        keys_by_type[result_type] = []
                    keys_by_type[result_type].append(f"{storage_path}.{key}")
        
        return keys_by_type
    
    @classmethod
    def _get_all_keys(cls, adata: AnnData) -> List[str]:
        """获取所有键"""
        all_keys = []
        
        for mod_name, storage_path in cls.MODULE_STORAGE_MAP.items():
            storage_data = cls._get_nested_data(adata.uns, storage_path)
            if storage_data:
                for key in storage_data.keys():
                    all_keys.append(key)
        
        return all_keys


# === 便捷函数 ===

def check_storage_status(
    adata: AnnData,
    module: Optional[str] = None,
    include_size: bool = True,
    include_details: bool = False
) -> Dict:
    """检查存储状态（便捷函数）"""
    return StorageManager.get_storage_status(
        adata, module, include_size, include_details
    )

def cleanup_storage(
    adata: AnnData,
    module: Optional[str] = None,
    result_types: Optional[List[str]] = None,
    keys_to_keep: Optional[List[str]] = None,
    dry_run: bool = False,
    force: bool = False,
    create_backup: bool = True,
    backup_dir: Optional[str] = None
) -> Dict:
    """
        Cleanup storage (convenience function)
        Parameters:
            create_backup: If True, save deleted data to backup file
            backup_dir: Directory for backups (default: ./sclucid_backups)
    """
    if create_backup and not dry_run:
        if backup_dir is None:
            backup_dir = Path.cwd() / "sclucid_backups"
        backup_dir = Path(backup_dir)
        backup_dir.mkdir(exist_ok=True)
        
        # 收集要删除的数据
        backup_data = {}
        for key_path in keys_to_delete:
            value = cls._get_nested_data(adata.uns, key_path)
            backup_data[key_path] = value
        
        # 保存备份
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_file = backup_dir / f"cleanup_backup_{timestamp}.pkl"
        
        with open(backup_file, 'wb') as f:
            pickle.dump({
                'timestamp': timestamp,
                'module': module,
                'result_types': result_types,
                'data': backup_data
            }, f)
        
        log.info(f"Backup saved to: {backup_file}")
        cleanup_report['backup_file'] = str(backup_file)
    return StorageManager.cleanup_storage(
        adata, module, result_types, keys_to_keep, dry_run, force
    )

def optimize_storage(
    adata: AnnData,
    strategy: str = 'integrated_only',
    keep_configs: bool = True,
    dry_run: bool = False
) -> Dict:
    """优化存储策略（便捷函数）"""
    return StorageManager.optimize_storage(
        adata, strategy, keep_configs, dry_run
    )