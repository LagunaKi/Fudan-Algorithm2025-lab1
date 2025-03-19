import numpy as np
import matplotlib.pyplot as plt

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def build_hash_table(S, max_len, min_len=100):
    """
    构建哈希表，仅存储长度在 [min_len, max_len] 的子串及其反向互补。
    """
    hash_table = {}
    n = len(S)
    for i in range(n):
        for k in range(min_len, min(max_len, n - i) + 1):
            sub = S[i:i+k]
            hash_table[sub] = i
            hash_table[reverse_complement(sub)] = i
    return hash_table

def find_repeats_optimized(S_reference, S_query, max_len=float('inf'), min_len=40):
    n, m = len(S_reference), len(S_query)
    hash_table = build_hash_table(S_reference, max_len, min_len)
    dp = np.full((m + 1), float('inf'))
    dp[0] = 0
    backtrack = [None] * (m + 1)

    for j in range(1, m + 1):
        # 默认处理：逐个字符比对（视为插入/错配，此处不处理）
        dp[j] = dp[j-1] + 1
        backtrack[j] = (j-1, None, False)
        # 检查重复片段（仅处理长度 ≥ min_len）
        for k in range(min_len, min(max_len, j) + 1):
            sub = S_query[j - k:j]
            if sub in hash_table:
                if dp[j - k] + 1 < dp[j]:
                    dp[j] = dp[j - k] + 1
                    pos_ref = hash_table[sub]
                    is_reverse = (sub != S_reference[pos_ref:pos_ref + k])
                    backtrack[j] = (j - k, pos_ref, is_reverse)

    # 回溯路径
    repeats = []
    j = m
    while j > 0:
        prev_j, pos_ref, is_reverse = backtrack[j]
        if pos_ref is not None:
            k = j - prev_j
            repeats.append((pos_ref, k, 1, is_reverse))
        j = prev_j

    return repeats[::-1]

# 绘制 Dot Plot
def plot_dot(S_reference, S_query):
    # 初始化坐标列表
    x_pos, y_pos = [], []     # 正向匹配点
    x_neg, y_neg = [], []     # 反向互补匹配点
    count = 0
    total = len(S_reference) * len(S_query)
    # 收集匹配点坐标
    for i in range(len(S_reference)):
        for j in range(len(S_query)):
            count += 1
            if S_reference[i] == S_query[j]:
                x_pos.append(j)
                y_pos.append(i)
            elif S_reference[i] == reverse_complement(S_query[j]):
                x_neg.append(j)
                y_neg.append(i)
        print(f"准备数据... {count/total*100:.2f}%")

    # 一次性绘制所有点
    plt.figure(figsize=(20, 12))
    plt.scatter(x_pos, y_pos, c='purple', s=1, alpha=0.5, label='Strand 1')
    plt.scatter(x_neg, y_neg, c='red', s=1, alpha=0.5, label='Strand -1')
    plt.title("Dot Plot")
    plt.xlabel("Query")
    plt.ylabel("Reference")
    plt.legend(loc='upper right')
    plt.xlim(0, len(S_query))
    plt.ylim(0, len(S_reference))
    plt.show()

# 测试示例
S_reference = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"
S_query = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"
print("优化后结果（限制重复长度）：")
repeats = find_repeats_optimized(S_reference, S_query, max_len=len(S_reference), min_len=50)
for idx, (pos, length, count, is_reverse) in enumerate(repeats, 1):
    print(f"重复{idx}: 位置={pos}, 长度={length}, 重复次数={count}, 反向={is_reverse}")


# 绘制 Dot Plot
print("\n绘制 Dot Plot...")
plot_dot(S_reference, S_query)
