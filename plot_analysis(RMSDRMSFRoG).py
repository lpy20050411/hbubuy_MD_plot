import matplotlib.pyplot as plt
import numpy as np
import os

# 读取name.txt文件获取复合物名称
def get_complex_name():
    try:
        with open('name.txt', 'r') as f:
            return f.read().strip()
    except:
        return 'Complex'

# 设置字体为Arial
plt.rcParams['font.family'] = 'Arial'

# 定义函数读取xvg文件
def read_xvg(file_path):
    data = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    x = float(parts[0])
                    y = float(parts[1])
                    data.append([x, y])
                except:
                    continue
    return np.array(data)

# 定义颜色映射
color_map = {
    'C-alpha_C-alpha': 'blue',
    'backbone_backbone': 'green',
    'backbone_ligand': 'brown',
    'backbone_protein-H': 'orange',
    'backbone_sidechain-H': 'red',
    'ligand_ligand': 'magenta',
    'C-alpha': 'blue',
    'backbone': 'green',
    'protein-H': 'orange',
    'sidechain-H': 'red'
}

# 绘制RMSD图
def plot_rmsd(ax, complex_name):
    # 获取RMSD文件
    rmsd_files = [f for f in os.listdir('.') if 'RMSD' in f and f.endswith('.xvg')]
    
    # 蓝色系颜色列表
    blue_colors = ['blue', 'deepskyblue', 'dodgerblue', 'mediumblue', 'navy']
    
    for i, rmsd_file in enumerate(rmsd_files):
        # 提取文件名中的类型
        file_type = rmsd_file.replace('.xvg', '')
        # 使用蓝色系颜色
        color = blue_colors[i % len(blue_colors)]
        
        # 读取数据
        data = read_xvg(rmsd_file)
        # 转换时间单位从ps到ns
        time = data[:, 0] / 1000
        rmsd = data[:, 1]
        
        # 绘制曲线
        ax.plot(time, rmsd, label=file_type, color=color)
    
    # 设置标题和标签
    ax.set_title(f'{complex_name} RMSD', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=10)
    ax.set_ylabel('RMSD (nm)', fontsize=10)
    
    # 添加图例
    ax.legend(loc='upper right', fontsize=8)
    
    # 添加小写序号在左上角
    ax.text(0.05, 0.95, '(a)', transform=ax.transAxes, fontsize=12, fontweight='bold')

# 绘制RMSF图
def plot_rmsf(ax, complex_name):
    # 获取RMSF文件
    rmsf_files = [f for f in os.listdir('.') if 'RMSF' in f and f.endswith('.xvg')]
    
    # 红色系颜色列表
    red_colors = ['red', 'crimson', 'firebrick', 'indianred', 'tomato']
    
    for i, rmsf_file in enumerate(rmsf_files):
        # 提取文件名中的类型
        file_type = rmsf_file.replace('.xvg', '')
        # 使用红色系颜色
        color = red_colors[i % len(red_colors)]
        
        # 读取数据
        data = read_xvg(rmsf_file)
        atom = data[:, 0]
        rmsf = data[:, 1]
        
        # 绘制曲线
        ax.plot(atom, rmsf, label=file_type, color=color)
    
    # 设置标题和标签
    ax.set_title(f'{complex_name} RMSF', fontsize=12)
    ax.set_xlabel('Atom', fontsize=10)
    ax.set_ylabel('RMSF (nm)', fontsize=10)
    
    # 添加图例
    ax.legend(loc='upper right', fontsize=8)
    
    # 添加小写序号在左上角
    ax.text(0.05, 0.95, '(b)', transform=ax.transAxes, fontsize=12, fontweight='bold')

# 绘制RoG图
def plot_rog(ax, complex_name):
    # 获取RoG文件
    rog_files = [f for f in os.listdir('.') if 'Gyrate' in f and f.endswith('.xvg')]
    
    # 绿色系颜色列表
    green_colors = ['green', 'limegreen', 'forestgreen', 'seagreen', 'mediumseagreen']
    
    for i, rog_file in enumerate(rog_files):
        # 提取文件名中的类型
        file_type = rog_file.replace('.xvg', '')
        # 使用绿色系颜色
        color = green_colors[i % len(green_colors)]
        
        # 读取数据
        data = read_xvg(rog_file)
        # 转换时间单位从ps到ns
        time = data[:, 0] / 1000
        rog = data[:, 1]
        
        # 绘制曲线
        ax.plot(time, rog, label=file_type, color=color)
    
    # 设置标题和标签
    ax.set_title(f'{complex_name} RoG', fontsize=12)
    ax.set_xlabel('Time (ns)', fontsize=10)
    ax.set_ylabel('RoG (nm)', fontsize=10)
    
    # 添加图例
    ax.legend(loc='upper right', fontsize=8)
    
    # 添加小写序号在左上角
    ax.text(0.05, 0.95, '(c)', transform=ax.transAxes, fontsize=12, fontweight='bold')

# 组合三张图
def plot_combined(complex_name):
    # 创建一个2x2的子图布局，第三张图占据整行
    fig = plt.figure(figsize=(12, 8))
    
    # 第一张图：左上角
    ax1 = fig.add_subplot(2, 2, 1)
    plot_rmsd(ax1, complex_name)
    
    # 第二张图：右上角
    ax2 = fig.add_subplot(2, 2, 2)
    plot_rmsf(ax2, complex_name)
    
    # 第三张图：下方整行
    ax3 = fig.add_subplot(2, 1, 2)
    plot_rog(ax3, complex_name)
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片，提高DPI到600
    plt.savefig('combined_plot.png', dpi=600)
    plt.close()

# 执行绘图
if __name__ == '__main__':
    # 获取复合物名称
    complex_name = get_complex_name()
    
    # 生成单独的图片
    # 创建临时图形对象来绘制单独的图
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    plot_rmsd(ax, complex_name)
    plt.tight_layout()
    plt.savefig('RMSD_plot.png', dpi=300)
    plt.close()
    
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    plot_rmsf(ax, complex_name)
    plt.tight_layout()
    plt.savefig('RMSF_plot.png', dpi=300)
    plt.close()
    
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    plot_rog(ax, complex_name)
    plt.tight_layout()
    plt.savefig('RoG_plot.png', dpi=300)
    plt.close()
    
    # 生成组合图片
    plot_combined(complex_name)
    
    print(f"绘图完成，生成了单独的图片和组合图片，使用复合物名称：{complex_name}")
