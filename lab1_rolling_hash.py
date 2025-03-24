import numpy as np
import matplotlib.pyplot as plt

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

class RollingHash:
    def __init__(self, base=128, mod=10**18+3):
        self.base = base
        self.mod = mod
        self.powers = [1]  # base^0 = 1
    
    def precompute_powers(self, max_len):
        while len(self.powers) <= max_len:
            self.powers.append((self.powers[-1] * self.base) % self.mod)
    
    def get_hash(self, s):
        hash_val = 0
        for c in s:
            hash_val = (hash_val * self.base + ord(c)) % self.mod
        return hash_val
    
    def update_hash(self, prev_hash, prev_char, new_char, k):
        return ((prev_hash - ord(prev_char) * self.powers[k-1]) * self.base + ord(new_char)) % self.mod

def build_hash_table(S, max_len, min_len=100):
    rh = RollingHash()
    rh.precompute_powers(max_len)
    
    hash_table = {}
    n = len(S)
    
    for k in range(min_len, max_len+1):
        if k > n: break

        # 预计算初始哈希
        current_hash = rh.get_hash(S[:k])
        hash_table.setdefault(current_hash, []).append((0, S[:k]))
        
        # 滑动窗口更新哈希
        for i in range(1, n - k + 1):
            prev_char = S[i-1]
            new_char = S[i+k-1]
            current_hash = rh.update_hash(current_hash, prev_char, new_char, k)
            sub = S[i:i+k]
            hash_table.setdefault(current_hash, []).append((i, sub))
            
            # 存储反向互补
            rc_sub = reverse_complement(sub)
            rc_hash = rh.get_hash(rc_sub)
            hash_table.setdefault(rc_hash, []).append((i, rc_sub))
    
    return hash_table

def find_repeats_optimized(S_reference, S_query, max_len=float('inf'), min_len=40):
    n, m = len(S_reference), len(S_query)
    hash_table = build_hash_table(S_reference, max_len, min_len)
    dp = np.full((m + 1), float('inf'))
    dp[0] = 0
    backtrack = [None] * (m + 1)
    rh = RollingHash()

    for j in range(1, m + 1):
        dp[j] = dp[j-1] + 1
        backtrack[j] = (j-1, None, False)
        
        for k in range(min_len, min(max_len, j) + 1):
            sub = S_query[j-k:j]
            sub_hash = rh.get_hash(sub)
            
            if sub_hash in hash_table:
                # 处理哈希冲突
                for pos_ref, stored_sub in hash_table[sub_hash]:
                    if stored_sub == sub:  # 实际字符串比对
                        if dp[j - k] + 1 < dp[j]:
                            dp[j] = dp[j - k] + 1
                            is_reverse = (sub != S_reference[pos_ref:pos_ref + k])
                            backtrack[j] = (j - k, pos_ref, is_reverse)
                        break  # 找到第一个匹配即停止

    # 回溯逻辑保持不变
    repeats = []
    j = m
    while j > 0:
        prev_j, pos_ref, is_reverse = backtrack[j]
        if pos_ref is not None:
            k = j - prev_j
            repeats.append((pos_ref, k, 1, is_reverse))
        j = prev_j

    return repeats[::-1]

# 测试示例
S_reference = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"
S_query = "CTGCAACGTTCGTGGTTCATGTTTGAGCGATAGGCCGAAACTAACCGTGCATGCAACGTTAGTGGATCATTGTGGAACTATAGACTCAAACTAAGCGAGCTTGCAACGTTAGTGGACCCTTTTTGAGCTATAGACGAAAACGGACCGAGGCTGCAAGGTTAGTGGATCATTTTTCAGTTTTAGACACAAACAAACCGAGCCATCAACGTTAGTCGATCATTTTTGTGCTATTGACCATATCTCAGCGAGCCTGCAACGTGAGTGGATCATTCTTGAGCTCTGGACCAAATCTAACCGTGCCAGCAACGCTAGTGGATAATTTTGTTGCTATAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCCTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCTAGACCAACACTAATCGAGACTGCCTCGTTAGTGCATCATTTTTGCGCCATAGACCATAGCTAAGCGAGCGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGGCTCGCTTAGCTATGGTCTATGGCGCAAAAATGATGCACTAACGAGGCAGTCTCGATTAGTGTTGGTCTATAGCAACAAAATTATCCACTAGCGTTGCTGCTTACCATCGGACCTCCACGAATCTGAAAAGTTTTAATTTCCGAGCGATACTTACGACCGGACCTCCACGAATCAGAAAGGGTTCACTATCCGCTCGATACATACGATCGGACCTCCACGACTCTGTAAGGTTTCAAAATCCGCACGATAGTTACGACCGTACCTCTACGAATCTATAAGGTTTCAATTTCCGCTGGATCCTTACGATCGGACCTCCTCGAATCTGCAAGGTTTCAATATCCGCTCAATGGTTACGGACGGACCTCCACGCATCTTAAAGGTTAAAATAGGCGCTCGGTACTTACGATCGGACCTCTCCGAATCTCAAAGGTTTCAATATCCGCTTGATACTTACGATCGCAACACCACGGATCTGAAAGGTTTCAATATCCACTCTATA"
print("优化后结果（限制重复长度）：")
repeats = find_repeats_optimized(S_reference, S_query, max_len=len(S_reference)//2, min_len=50)
for idx, (pos, length, count, is_reverse) in enumerate(repeats, 1):
    print(f"重复{idx}: 位置={pos}, 长度={length}, 重复次数={count}, 反向={is_reverse}")

