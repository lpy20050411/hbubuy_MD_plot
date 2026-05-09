import MDAnalysis as mda
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
import numpy as np
import sys
import multiprocessing as mp
from functools import partial



'''
做相互作用分析

'''
# 尝试导入pandas，如果不可用则使用基本文件操作
try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("pandas not available, using basic file operations for CSV output")

# 定义单帧分析函数
def analyze_frame(frame_idx, config):
    """分析单个帧中的配体-蛋白相互作用"""
    # 在子进程中独立加载轨迹文件
    import MDAnalysis as mda
    u = mda.Universe(config['gro_file'], config['xtc_file'])
    ligand = u.select_atoms(config['ligand_selection'])
    protein = u.select_atoms(config['protein_selection'])
    
    # 移动到指定帧
    u.trajectory[frame_idx]
    
    frame_interactions = {}
    
    # 为当前帧计算配体与蛋白残基的最小距离
    for residue in protein.residues:
        res_name = f"{residue.resname}_{residue.resid}"
        
        # 选择残基的所有原子
        residue_atoms = residue.atoms
        # 计算残基与配体之间的最小距离
        min_distance = float('inf')
        for lig_atom in ligand:
            for res_atom in residue_atoms:
                distance = np.linalg.norm(lig_atom.position - res_atom.position)
                if distance < min_distance:
                    min_distance = distance
        
        # 只考虑距离小于cutoff的残基
        if min_distance < config['distance_cutoff']:
            # 初始化残基的相互作用计数器
            if res_name not in frame_interactions:
                frame_interactions[res_name] = {'HydrogenBond': 0, 'Hydrophobic': 0, 'VanDerWaals': 0}
            
            # 基于距离和残基类型分配相互作用类型
            if min_distance < config['hydrogen_bond_cutoff']:
                # 氢键或离子键
                if residue.resname in ['ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'SER', 'THR', 'TYR', 'ASN', 'GLN']:
                    frame_interactions[res_name]['HydrogenBond'] = 1
                else:
                    frame_interactions[res_name]['VanDerWaals'] = 1
            elif min_distance < config['hydrophobic_cutoff']:
                # 疏水相互作用
                if residue.resname in ['ALA', 'VAL', 'LEU', 'ILE', 'PHE', 'TRP', 'MET', 'PRO']:
                    frame_interactions[res_name]['Hydrophobic'] = 1
                else:
                    frame_interactions[res_name]['VanDerWaals'] = 1
            else:
                # 范德华相互作用
                frame_interactions[res_name]['VanDerWaals'] = 1
    
    return frame_interactions

if __name__ == '__main__':
    # 用户配置区
    CONFIG = {
        # 输入文件
        'gro_file': 'md.gro',
        'xtc_file': 'md.xtc',
        
        # 配体选择
        'ligand_selection': 'resname MOL',
        'protein_selection': 'protein',
        
        # 相互作用参数
        'distance_cutoff': 4.0,  # 残基与配体的最大距离
        'hydrogen_bond_cutoff': 2.5,  # 氢键距离阈值
        'hydrophobic_cutoff': 3.5,  # 疏水相互作用距离阈值
        
        # 多线程设置
        'use_multithreading': True,
        'num_threads': min(mp.cpu_count(), 64),  # 限制线程数不超过64
        
        # 输出设置
        'output_plot': 'interaction_frequencies.png',
        'output_csv': 'interaction_frequencies.csv',
        'dpi': 300,
        'figsize': (15, 6),
        'font_size': 12,
        'xtick_rotation': 45
    }

    print("Starting analysis...")
    print(f"Using configuration: {CONFIG}")

    # 加载轨迹文件
    try:
        print("Loading trajectory files...")
        u = mda.Universe(CONFIG['gro_file'], CONFIG['xtc_file'])
        print(f"Trajectory loaded successfully. Number of frames: {len(u.trajectory)}")
    except Exception as e:
        print(f"Error loading trajectory: {e}")
        sys.exit(1)

    # 选择配体和蛋白
    try:
        print("Selecting ligand and protein...")
        ligand = u.select_atoms(CONFIG['ligand_selection'])
        protein = u.select_atoms(CONFIG['protein_selection'])
        print(f"Ligand atoms: {len(ligand)}")
        print(f"Protein atoms: {len(protein)}")
    except Exception as e:
        print(f"Error selecting atoms: {e}")
        sys.exit(1)

    # 计算配体与蛋白残基的距离
    print("Calculating distances between ligand and protein residues...")

    # 初始化相互作用计数器
    interactions = {}
    total_frames = len(u.trajectory)

    if CONFIG['use_multithreading']:
        print(f"Using multiprocessing with {CONFIG['num_threads']} threads...")
        
        # 创建部分函数，只传递配置参数
        analyze_frame_partial = partial(analyze_frame, config=CONFIG)
        
        # 使用多进程池处理所有帧
        with mp.Pool(processes=CONFIG['num_threads']) as pool:
            # 处理所有帧
            frame_results = pool.map(analyze_frame_partial, range(total_frames))
            
            # 合并所有帧的结果
            for i, frame_interactions in enumerate(frame_results):
                if i % 10 == 0:
                    print(f"Processing frame {i}/{total_frames}")
                
                for res_name, res_interactions in frame_interactions.items():
                    if res_name not in interactions:
                        interactions[res_name] = {'HydrogenBond': 0, 'Hydrophobic': 0, 'VanDerWaals': 0}
                    
                    for interaction_type, count in res_interactions.items():
                        interactions[res_name][interaction_type] += count
    else:
        print("Using single-threaded processing...")
        # 单线程处理
        for i, ts in enumerate(u.trajectory):
            if i % 10 == 0:
                print(f"Processing frame {i}/{total_frames}")
            
            frame_interactions = analyze_frame(i, CONFIG)
            
            for res_name, res_interactions in frame_interactions.items():
                if res_name not in interactions:
                    interactions[res_name] = {'HydrogenBond': 0, 'Hydrophobic': 0, 'VanDerWaals': 0}
                
                for interaction_type, count in res_interactions.items():
                    interactions[res_name][interaction_type] += count

    print(f"Found {len(interactions)} residues with interactions across {total_frames} frames")

    # 转换为频率
    print("Calculating interaction frequencies...")
    interaction_freq = {}
    for res, counts in interactions.items():
        interaction_freq[res] = {}
        for interaction_type, count in counts.items():
            interaction_freq[res][interaction_type] = count / total_frames

    # 准备数据用于绘图
    print("Preparing data for plotting...")
    residues = list(interaction_freq.keys())
    interaction_types = set()
    freq_data = {}

    for res, types in interaction_freq.items():
        for interaction_type in types:
            interaction_types.add(interaction_type)
            if interaction_type not in freq_data:
                freq_data[interaction_type] = []
            freq_data[interaction_type].append(types.get(interaction_type, 0))

    # 确保所有残基都有所有相互作用类型的数据
    for interaction_type in interaction_types:
        while len(freq_data[interaction_type]) < len(residues):
            freq_data[interaction_type].append(0)

    # 排序残基
    residues.sort()
    print(f"Residues: {residues}")
    print(f"Interaction types: {interaction_types}")

    # 定义颜色映射
    color_map = {
        'HydrogenBond': '#4CAF50',  # 绿色
        'Hydrophobic': '#2196F3',    # 蓝色
        'PiStacking': '#FF9800',      # 橙色
        'Ionic': '#F44336',           # 红色
        'VanDerWaals': '#9C27B0'      # 紫色
    }

    # 绘制堆叠柱状图
    print("Plotting results...")
    try:
        plt.figure(figsize=CONFIG['figsize'])
        plt.rcParams.update({'font.size': CONFIG['font_size']})
        
        bottom = np.zeros(len(residues))
        for interaction_type in interaction_types:
            if interaction_type in color_map:
                color = color_map[interaction_type]
            else:
                color = '#607D8B'  # 默认颜色
            
            plt.bar(residues, freq_data[interaction_type], bottom=bottom, label=interaction_type, color=color)
            bottom += np.array(freq_data[interaction_type])

        plt.xlabel('Residue')
        plt.ylabel('Interaction Fraction')
        plt.title('EGFR-Lutolin Complex Interaction Frequencies')
        plt.xticks(rotation=CONFIG['xtick_rotation'], ha='right')
        plt.legend()
        plt.tight_layout()
        plt.savefig(CONFIG['output_plot'], dpi=CONFIG['dpi'])
        print(f"Plot saved as {CONFIG['output_plot']}")
    except Exception as e:
        print(f"Error plotting: {e}")
        import traceback
        traceback.print_exc()

    # 保存数据到CSV
    try:
        if HAS_PANDAS:
            # 使用pandas写入CSV
            df = pd.DataFrame(freq_data, index=residues)
            df.to_csv(CONFIG['output_csv'])
            print(f"Data saved as {CONFIG['output_csv']} (using pandas)")
        else:
            # 使用基本的Python写入CSV
            with open(CONFIG['output_csv'], 'w') as f:
                # 写入表头
                header = ['Residue'] + list(interaction_types)
                f.write(','.join(header) + '\n')
                
                # 写入数据
                for i, res in enumerate(residues):
                    row = [res]
                    for interaction_type in interaction_types:
                        row.append(str(freq_data[interaction_type][i]))
                    f.write(','.join(row) + '\n')
            print(f"Data saved as {CONFIG['output_csv']} (using basic file operations)")
    except Exception as e:
        print(f"Error saving CSV: {e}")
        import traceback
        traceback.print_exc()

    print("Analysis completed!")
