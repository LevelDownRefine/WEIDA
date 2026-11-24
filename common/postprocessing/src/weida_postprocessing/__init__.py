"""共用后处理的统计入口。"""
from .author_statistics import (
    CONTRASTS,
    GROUPS,
    compare_groups,
    describe_groups,
    get_fdr,
    group_comparison,
)

__all__ = ['CONTRASTS', 'GROUPS', 'compare_groups', 'describe_groups', 'get_fdr', 'group_comparison']
