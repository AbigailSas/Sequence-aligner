import numpy as np

Bases = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
DIAG = 1
UP = 2
LEFT = 3


class Overlap:

    def __init__(self, score_matrix, seq_a, seq_b):
        self.score_matrix = score_matrix
        self.seq_a = seq_a
        self.seq_b = seq_b
        self.align_mat = np.zeros((len(seq_a) + 1, len(seq_b) + 1))
        self.pointer_mat = np.zeros((len(seq_a) + 1, len(seq_b) + 1))
        self.score = 0

    def get_score(self, a, b):
        return self.score_matrix[Bases[a], Bases[b]]

    def init_overlap(self):
        """
        initialize mat(i,0), mat(0,j) to be 0
        :return: void function
        """
        sum_of_gaps = 0
        for i in range(1, len(self.seq_a) + 1):  # init pointer matrix
            self.pointer_mat[i, 0] = LEFT
        for i in range(1, len(self.seq_b) + 1):
            score = self.get_score("-", self.seq_b[i - 1])
            self.align_mat[0, i] = score + sum_of_gaps
            sum_of_gaps += score
            self.pointer_mat[0, i] = UP

    def fill_overlap(self):
        self.init_overlap()  # initialize matrix

        for j in range(1, len(self.seq_b) + 1):  # fill matrix
            for i in range(1, len(self.seq_a) + 1):
                diag = self.align_mat[i - 1, j - 1] + self.get_score(self.seq_a[i - 1], self.seq_b[j - 1])
                left = self.align_mat[i - 1, j] + self.get_score(self.seq_a[i - 1], "-")
                up = self.align_mat[i, j - 1] + self.get_score("-", self.seq_b[j - 1])
                self.align_mat[i, j] = max(diag, up, left)
                self.align_mat[i, j] = max(diag, up, left)
                if max(diag, up, left) == diag:
                    self.pointer_mat[i, j] = DIAG
                elif max(diag, up, left) == up:
                    self.pointer_mat[i, j] = UP
                else:
                    self.pointer_mat[i, j] = LEFT
        self.score = self.align_mat[len(self.seq_a), len(self.seq_b)]

    def overlap_get_solution(self):
        final_a, final_b = "", ""
        i = len(self.seq_a)
        j = np.argmax(self.align_mat[-1, :])
        j_final = j
        while i > 0 or j > 0:
            if self.pointer_mat[i, j] == DIAG:
                final_a = self.seq_a[i-1] + final_a
                final_b = self.seq_b[j-1] + final_b
                i -= 1
                j -= 1
            elif self.pointer_mat[i, j] == LEFT:  # mat[i - 1, j] + score(seq_a[i], "-")
                final_a = self.seq_a[i-1] + final_a
                final_b = "-" + final_b
                i -= 1
            else:  # UP
                final_a = "-" + final_a
                final_b = self.seq_b[j-1] + final_b
                j -= 1

        if j_final != len(self.seq_b):
            while j_final < len(self.seq_b):
                final_a = final_a + "-"
                final_b = final_b + self.seq_b[j_final]
                j_final += 1
        return final_a, final_b

    def get_alignment(self):
        self.fill_overlap()
        return self.overlap_get_solution()

    def get_align_score(self):
        return self.score
