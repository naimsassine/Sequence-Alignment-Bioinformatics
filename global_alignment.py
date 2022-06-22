# parse the substitution matrices 
# read from file into a dict
import pdb
def subs_parsing_blosum(file_name) :
    matrix = []
    l = []
    nonos = ["#", "", " ", "*"]
    with open(file_name, 'r') as f:
        for line in f : 
            if line[0] not in nonos and line != "\n": 
                l.append(line[2:71].split(' ')) 

    for element in l : 
        element = [int(i) for i in element if i != '' and i != '\n']
        matrix.append(element)
    
    return matrix


def subs_parsing_pam(file_name) :
    matrix = []
    l = []
    nonos = ["#", "", " ", "*"]
    with open(file_name, 'r') as f:
        for line in f : 
            if line[0] not in nonos and line != "\n": 
                l.append(line[2:93].split(' ')) 

    for element in l : 
        element = [int(i) for i in element if i != '' and i != '\n']
        matrix.append(element)
    
    return matrix

# A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X

def matrix_to_dict(matrix) :
    aa_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X"] 
    dico = {}
    for i in range(0, len(matrix)) : 
        for j in range(0, len(matrix[0])) : 
            dico[aa_list[i] + aa_list[j]] = matrix[i][j]
    return dico



# get the dictionnaries containing the substitution matrices
dico_blosum62 = matrix_to_dict(subs_parsing_blosum("data/blosum62.txt"))
dico_blosum80 = matrix_to_dict(subs_parsing_blosum("data/blosum80.txt"))
dico_pam60 = matrix_to_dict(subs_parsing_pam("data/pam60.txt"))
dico_pam120 = matrix_to_dict(subs_parsing_blosum("data/pam120.txt"))


# Needleman-Wunsch algorithm
# how to find the K best solutions ? 
# f there are multiple arrows to choose from, they represent a branching of the alignments. 
# If two or more branches all belong to paths from the bottom right to the top left cell, 
# they are equally viable alignments. In this case, note the paths as separate alignment candidates.


class SequenceAlignment(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.solution = []
        self.aligned_seq = [[],[]]
        
    delta = lambda self, x, y, i, j: 1 if x[i] != y[j] else 0

    def locate_min(a):
        smallest = min(a)
        return smallest, [index for index, element in enumerate(a) 
                          if smallest == element]

    def find_solution(self, OPT, m, n, solution = [], aligned_seq = [[],[]]):

        # pdb.set_trace()
        if m == 0 and n == 0:
            print(solution, aligned_seq) 
            return
        else : 
            # We can only do insert if n != 0, align if there are element in both x, y, etc.
            insert = OPT[m][n - 1] + 1 if n != 0 else float("inf")
            align = (
                OPT[m - 1][n - 1] + self.delta(self.x, self.y, m - 1, n - 1)
                if m != 0 and n != 0
                else float("inf")
            )
            delete = OPT[m - 1][n] + 1 if m != 0 else float("inf")

            best_choice = min(insert, align, delete)

            if best_choice == insert and best_choice == align : 
                solution.append("insert_" + str(self.y[n - 1]))
                aligned_seq[1].insert(0, str(self.y[n - 1]))
                aligned_seq[0].insert(0, "-")
                self.find_solution(OPT, m, n - 1, solution, aligned_seq)
                # second branch
                sb_solution = solution.copy()
                sb_aligned_seq = aligned_seq.copy()
                sb_solution.append("align_" + str(self.y[n - 1]))
                sb_aligned_seq[1].insert(0, str(self.y[n - 1]))
                sb_aligned_seq[0].insert(0, str(self.y[n - 1]))
                self.find_solution(OPT, m - 1, n - 1, sb_solution, sb_aligned_seq)
            
            elif best_choice == insert and best_choice == delete: 
                solution.append("insert_" + str(self.y[n - 1]))
                aligned_seq[1].insert(0, str(self.y[n - 1]))
                aligned_seq[0].insert(0, "-")
                self.find_solution(OPT, m, n - 1, solution, aligned_seq)
                # second branch
                sb_solution = solution.copy()
                sb_aligned_seq = aligned_seq.copy()
                sb_solution.append("remove_" + str(self.x[m - 1]))
                sb_aligned_seq[0].insert(0, str(self.x[m - 1]))
                sb_aligned_seq[1].insert(0, "-")
                self.find_solution(OPT, m - 1, n - 1, sb_solution, sb_aligned_seq)
                
            elif best_choice == align and best_choice == delete : 
                solution.append("align_" + str(self.y[n - 1]))
                aligned_seq[1].insert(0, str(self.y[n - 1]))
                aligned_seq[0].insert(0, str(self.y[n - 1]))
                self.find_solution(OPT, m - 1, n - 1, solution, aligned_seq)
                # second branch
                sb_solution = solution.copy()
                sb_aligned_seq = aligned_seq.copy()
                sb_solution.append("remove_" + str(self.x[m - 1]))
                sb_aligned_seq[0].insert(0, str(self.x[m - 1]))
                sb_aligned_seq[1].insert(0, "-")
                self.find_solution(OPT, m - 1, n - 1, sb_solution, sb_aligned_seq)
                
            elif best_choice == align and best_choice == delete and best_choice == insert : 
                solution.append("insert_" + str(self.y[n - 1]))
                aligned_seq[1].insert(0, str(self.y[n - 1]))
                aligned_seq[0].insert(0, "-")
                self.find_solution(OPT, m, n - 1, solution, aligned_seq)
                # second branch
                sb_solution = solution.copy()
                sb_aligned_seq = aligned_seq.copy()
                sb_solution.append("align_" + str(self.y[n - 1]))
                sb_aligned_seq[1].insert(0, str(self.y[n - 1]))
                sb_aligned_seq[0].insert(0, str(self.y[n - 1]))
                self.find_solution(OPT, m - 1, n - 1, sb_solution, sb_aligned_seq)
                sb2_solution = solution.copy()
                sb2_aligned_seq = aligned_seq.copy()
                sb2_solution.append("remove_" + str(self.x[m - 1]))
                sb2_aligned_seq[0].insert(0, str(self.x[m - 1]))
                sb2_aligned_seq[1].insert(0, "-")
                self.find_solution(OPT, m - 1, n - 1, sb2_solution, sb2_aligned_seq)
                
            elif best_choice == insert :
                solution.append("insert_" + str(self.y[n - 1]))
                aligned_seq[1].insert(0, str(self.y[n - 1]))
                aligned_seq[0].insert(0, "-")
                self.find_solution(OPT, m, n - 1, solution, aligned_seq)

            elif best_choice == align:
                solution.append("align_" + str(self.y[n - 1]))
                aligned_seq[1].insert(0, str(self.y[n - 1]))
                aligned_seq[0].insert(0, str(self.y[n - 1]))
                self.find_solution(OPT, m - 1, n - 1, solution, aligned_seq)

            elif best_choice == delete:
                solution.append("remove_" + str(self.x[m - 1]))
                aligned_seq[0].insert(0, str(self.x[m - 1]))
                aligned_seq[1].insert(0, "-")
                self.find_solution(OPT, m - 1, n, solution, aligned_seq)
            

            
        

        
    def alignment(self):
        n = len(self.y)
        m = len(self.x)
        OPT = [[0 for i in range(n + 1)] for j in range(m + 1)]

        for i in range(1, m + 1):
            OPT[i][0] = i

        for j in range(1, n + 1):
            OPT[0][j] = j

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                OPT[i][j] = min(
                    OPT[i - 1][j - 1] + self.delta(self.x, self.y, i - 1, j - 1),
                    OPT[i - 1][j] + 1,
                    OPT[i][j - 1] + 1,
                )  # align, delete, insert respectively

        self.find_solution(OPT, m, n)


    
x = 'TGACGTGC'
y = 'TCGACGTCA'
print('We we want to transform: ' + x + ' to: ' + y)
sqalign = SequenceAlignment(x, y)
sqalign.alignment()

# how to get the final aligment : 