from collections import defaultdict

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

class State:
    def __init__(self):
        self.length = 0
        self.link = -1
        self.transitions = {}

class SAM:
    def __init__(self):
        self.size = 1
        self.last = 0
        self.states = [State()]
    
    def sa_extend(self, c):
        p = self.last
        curr = self.size
        self.size += 1
        self.states.append(State())
        self.states[curr].length = self.states[p].length + 1
        while p >= 0 and c not in self.states[p].transitions:
            self.states[p].transitions[c] = curr
            p = self.states[p].link
        if p == -1:
            self.states[curr].link = 0
        else:
            q = self.states[p].transitions[c]
            if self.states[p].length + 1 == self.states[q].length:
                self.states[curr].link = q
            else:
                clone = self.size
                self.size += 1
                self.states.append(State())
                self.states[clone].length = self.states[p].length + 1
                self.states[clone].transitions = self.states[q].transitions.copy()
                self.states[clone].link = self.states[q].link
                while p >= 0 and self.states[p].transitions[c] == q:
                    self.states[p].transitions[c] = clone
                    p = self.states[p].link
                self.states[q].link = clone
                self.states[curr].link = clone
        self.last = curr

def build_sam(s):
    sam = SAM()
    for c in s:
        sam.sa_extend(c)
    return sam

def compute_match_lengths(s, sam, reverse=False):
    n = len(s)
    current_state = 0
    current_len = 0
    match_lengths = [0] * (n + 1)
    for i in range(n):
        c = s[i]
        if reverse:
            c = reverse_complement(c)
        while current_state != -1 and c not in sam.states[current_state].transitions:
            current_state = sam.states[current_state].link
            current_len = sam.states[current_state].length if current_state != -1 else 0
        if current_state == -1:
            current_state = 0
            current_len = 0
        else:
            current_state = sam.states[current_state].transitions[c]
            current_len += 1
        match_lengths[i+1] = current_len
    return match_lengths

def find_repeats(s_ref, s_query, min_len=50):
    s_ref_rev = reverse_complement(s_ref)
    sam_forward = build_sam(s_ref)
    sam_reverse = build_sam(s_ref_rev)
    
    forward_match = compute_match_lengths(s_query, sam_forward)
    reverse_match = compute_match_lengths(s_query, sam_reverse)
    
    m = len(s_query)
    dp = [float('inf')] * (m + 1)
    dp[0] = 0
    backtrack = [None] * (m + 1)
    
    for j in range(1, m + 1):
        dp[j] = dp[j-1] + 1
        backtrack[j] = (j-1, None, None)
        
        k = forward_match[j]
        if k >= min_len:
            if dp[j - k] + 1 < dp[j]:
                dp[j] = dp[j - k] + 1
                backtrack[j] = (j - k, k, False)
        
        k = reverse_match[j]
        if k >= min_len:
            if dp[j - k] + 1 < dp[j]:
                dp[j] = dp[j - k] + 1
                backtrack[j] = (j - k, k, True)
    
    repeats = []
    j = m
    while j > 0:
        prev_j, k, is_reverse = backtrack[j]
        if k is not None:
            start = prev_j
            repeat_seq = s_query[start:start + k]
            # 查找在reference中的位置
            if is_reverse:
                target_seq = reverse_complement(repeat_seq)
            else:
                target_seq = repeat_seq
            pos = s_ref.find(target_seq)
            if pos == -1:
                j = prev_j
                continue  # 未找到，跳过
            
            # 计算重复次数
            count = 1
            current = start + k
            while current + k <= m and s_query[current:current + k] == repeat_seq:
                count += 1
                current += k
            
            repeats.append((pos, k, count - 1, is_reverse))
        j = prev_j
    
    return repeats[::-1]

# 读取 reference 和 query 序列的函数
def read_reference(file_path):
    with open(file_path, 'r') as file:
        return file.read().replace('\n', '')

def read_query(file_path):
    with open(file_path, 'r') as file:
        return file.read().replace('\n', '')

# 测试示例
S_reference = read_reference('reference.txt')
S_query = read_query('query.txt')
repeats = find_repeats(S_reference, S_query, min_len=50)
for idx, (pos, length, count, is_reverse) in enumerate(repeats, 1):
    print(f"重复{idx}: 位置={pos}, 长度={length}, 重复次数={count}, 反向={is_reverse}")