# 导入必要的库
import numpy as np
import matplotlib.pyplot as plt

# 定义反向互补函数
def reverse_complement(seq):
    """
    计算 DNA 序列的反向互补。
    输入: DNA 序列 (例如 "ATCG")
    输出: 反向互补序列 (例如 "CGAT")
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# 初始化 DP 表格
def initialize_dp(n, m):
    """
    初始化 DP 表格，维度为 (n+1) x (m+1)。
    输入: n (S_reference 的长度), m (S_query 的长度)
    输出: 初始化后的 DP 表格
    """
    dp = np.full((n+1, m+1), float('inf'))  # 用无穷大初始化
    dp[0][0] = 0  # 起点权重为 0
    return dp

# 路径回溯函数
def backtrack_path(dp, S_reference, S_query):
    """
    从 DP 表格回溯路径，识别重复片段。
    输入: DP 表格, 参考序列 S_reference, 查询序列 S_query
    输出: 重复片段信息列表 [(位置, 长度, 重复次数, 是否反向), ...]
    """
    path = []
    i, j = len(S_reference), len(S_query)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and dp[i][j] == dp[i-1][j-1] and S_reference[i-1] == S_query[j-1]:
            # 如果是通过匹配到达的，继续回溯
            i -= 1
            j -= 1
        else:
            # 检查是否存在重复片段
            for k in range(1, j+1):
                if dp[i][j] == dp[i][j-k] + 1:
                    # 检查正向重复或反向重复
                    for p in range(len(S_reference) - k + 1):
                        if S_query[j-k:j] == S_reference[p:p+k]:
                            path.append((p, k, 1, False))  # 正向重复
                            j -= k
                            break
                        elif S_query[j-k:j] == reverse_complement(S_reference[p:p+k]):
                            path.append((p, k, 1, True))  # 反向重复
                            j -= k
                            break
                    else:
                        continue
                    break
            else:
                # 如果没有找到重复，可能是其他操作（本实验未处理）
                break
    return path[::-1]  # 反转路径以按顺序输出
# 绘制 Dot Plot
def plot_dot(S_reference, S_query):
    plt.figure(figsize=(8, 4))
    plt.title("Dot Plot")
    plt.xlabel("Query")
    plt.ylabel("Reference")
    plt.grid(True)
    
    # 找到匹配点
    for i in range(len(S_reference)):
        for j in range(len(S_query)):
            if S_reference[i] == S_query[j]:
                plt.scatter(j, i, c='purple', s=20, label='Strand 1' if i == 0 and j == 0 else "")
            elif S_reference[i] == reverse_complement(S_query[j]):
                plt.scatter(j, i, c='red', s=20, label='Strand -1' if i == 0 and j == 0 else "")
    
    plt.legend(loc='upper right')
    plt.xlim(0, len(S_query))
    plt.ylim(0, len(S_reference))
    plt.show()
# 主函数：识别重复片段并输出结果
def find_repeats(S_reference, S_query):
    """
    使用动态规划识别 S_query 中相对于 S_reference 的重复片段。
    输入: 参考序列 S_reference, 查询序列 S_query
    输出: 打印重复片段的位置、长度、重复次数和是否为反向重复
    """
    n, m = len(S_reference), len(S_query)
    dp = initialize_dp(n, m)
    
    # 填充 DP 表格
    for i in range(n+1):
        for j in range(m+1):
            if i > 0 and j > 0 and S_reference[i-1] == S_query[j-1]:
                # 如果当前字符匹配，则权重不变
                dp[i][j] = min(dp[i][j], dp[i-1][j-1])
            # 检查是否存在重复片段
            for k in range(1, j+1):
                for p in range(n - k + 1):
                    if (S_query[j-k:j] == S_reference[p:p+k] or 
                        S_query[j-k:j] == reverse_complement(S_reference[p:p+k])):
                        dp[i][j] = min(dp[i][j], dp[i][j-k] + 1)
    
    # 回溯路径并获取重复片段
    repeats = backtrack_path(dp, S_reference, S_query)
    
    # 输出结果
    print(f"参考序列: {S_reference}")
    print(f"查询序列: {S_query}")
    print("重复片段信息:")
    for repeat in repeats:
        pos, length, count, is_reverse = repeat
        print(f"位置: {pos}, 长度: {length}, 重复次数: {count}, 是否反向: {is_reverse}")
    
    # 绘制 Dot Plot
    plot_dot(S_reference, S_query)

# 测试示例
S_reference = "ATCG"
S_query = "ATCGTAGC"
print("测试示例：")
find_repeats(S_reference, S_query)