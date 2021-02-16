import hierarchicalCluster
import random
import timeit


def euclidean_distance(p1, p2):
    return hierarchicalCluster.euclidean_distance(p1, p2)


def manhattan_distance(p1, p2):
    return hierarchicalCluster.manhattan_distance(p1, p2)


def read_gene_expression_profiles(filename):
    rows = [line.rstrip().split("\t") for line in open(filename)]
    sample_names = rows[0]
    columns = zip(*rows[1:])
    profiles = [tuple(map(float, column)) for column in columns]
    return profiles, sample_names


cell_type_profiles, cell_type_names = read_gene_expression_profiles("cell_type_expression.txt")
blood_cell_type_names = [line.rstrip() for line in open("blood_cell_types.txt")]


# we will generate a random dataset of 100-dimensional profiles from the vertices of a hypercube
def random_hypercube_vertex(dims):
    return tuple(random.randint(0, 1) for i in range(dims))


random.seed(42)
dims = 100
num_profiles = 250
random_100_250_profiles = [random_hypercube_vertex(dims) for i in range(num_profiles)]
random_100_250_names = ["P{}".format(i) for i in range(num_profiles)]

# a dictionary of datasets for testing
# the key is the name of the test
# the value is a tuple (profiles, profile_names)
datasets = {
    "pair": ([(4, 2), (-2, -6)],
             ["A", "B"]),
    "triple": ([(4, 2), (-2, -6), (1, 6)],
               ["A", "B", "C"]),
    "quintet": ([(0,), (6,), (8,), (11,), (15,)],
                ["A", "B", "C", "D", "E"]),
    "tiebreaker": ([(0,), (1,), (2,), (3,), (4,), (5,)],
                   ["A", "B", "C", "D", "E", "F"]),
    "cell_type_sample": (
        [(2.7, 0.0, 2.2, 0.0, 1.6, 2.2, 0.8, 2.6, 0.0, 0.0),
         (3.0, 0.0, 2.0, 0.0, 3.2, 1.6, 1.2, 2.6, 1.1, 1.1),
         (2.9, 0.0, 2.3, 0.0, 2.2, 1.9, 0.5, 2.3, 1.0, 0.3),
         (2.6, 0.0, 2.3, 0.0, 2.2, 2.0, 0.7, 2.0, 1.3, 0.8),
         (2.7, 0.0, 2.1, 0.0, 4.0, 2.1, 1.3, 2.6, 1.0, 1.0),
         (2.8, 0.0, 2.1, 0.0, 2.6, 2.2, 2.3, 2.4, 1.4, 1.1),
         (2.7, 0.0, 2.1, 0.0, 3.2, 2.8, 1.9, 1.9, 1.2, 0.8),
         (2.8, 0.0, 2.5, 1.3, 2.8, 2.1, 2.8, 2.0, 1.0, 0.7),
         (2.9, 0.0, 2.0, 0.8, 1.4, 2.2, 1.9, 2.6, 1.0, 0.7)],
        ["placental pericyte",
         "stromal cell",
         "pericyte cell",
         "skin fibroblast",
         "hematopoietic cell",
         "stromal cell of ovary",
         "calvarial osteoblast",
         "osteoblast",
         "astrocyte"]),
    "random_100_250": (random_100_250_profiles, random_100_250_names),
    "cell_type": (cell_type_profiles, cell_type_names)
}


# testing functions
def test_case_newick(name, linkage, dist):
    profiles, names = datasets[name]
    tree = hierarchicalCluster.cluster_bottom_up(profiles, names, linkage=linkage, distance=dist)
    tree.sort_descendants()
    return tree.write(format=1)


def test_case(name, linkage, distance, correct_newick):
    output_newick = test_case_newick(name, linkage, distance)
    if output_newick != correct_newick:
        assert False, "Failed test\n Output: %s\nCorrect: %s" % (output_newick, correct_newick)
    else:
        print("SUCCESS: test passed")


# test pair_single_euclidean
test_case("pair", "single", euclidean_distance, "(A:5,B:5);")

# test pair_single_manhattan
test_case("pair", "single", manhattan_distance, "(A:7,B:7);")

# test triple_single_manhattan
test_case("triple", "single", manhattan_distance, "((A:3.5,C:3.5):3.5,B:7);")

# test triple_complete_manhattan
test_case("triple", "complete", manhattan_distance, "((A:3.5,C:3.5):4,B:7.5);")

# test triple_average_manhattan
test_case("triple", "average", manhattan_distance, "((A:3.5,C:3.5):3.75,B:7.25);")

# test quintet_single_manhattan
test_case("quintet", "single", manhattan_distance, "(A:3,(((B:1,C:1):0.5,D:1.5):0.5,E:2):1);")

# test quintet_complete_manhattan
test_case("quintet", "complete", manhattan_distance, "((A:4,(B:1,C:1):3):3.5,(D:2,E:2):5.5);")

# test quintet_average_manhattan
test_case("quintet", "average", manhattan_distance, "(A:5,((B:1,C:1):2,(D:2,E:2):1):2);")

# test tiebreaker_single_manhattan
test_case("tiebreaker", "single", manhattan_distance,
          "(((A:0.5,B:0.5):0,(C:0.5,D:0.5):0):0,(E:0.5,F:0.5):0);")

# test cell_type_sample_average_euclidean
test_case("cell_type_sample", "average", euclidean_distance,
          "((astrocyte:0.936623,((pericyte cell:0.377492,skin fibroblast:0.377492):0.3967,"
          "placental pericyte:0.774191):0.162432):0.172292,(((calvarial osteoblast:0.563471,"
          "stromal cell of ovary:0.563471):0.223098,"
          "(hematopoietic cell:0.504975,stromal cell:0.504975):0.281594):0.24476,osteoblast:1.03133):0.0775861);")

# test cell_type_sample_single_chebyshev
test_case("cell_type_sample",
          "single",
          lambda p1, p2: max(abs(e1 - e2) for e1, e2 in zip(p1, p2)),  # Chebyshev distance
          '((astrocyte:0.55,(((calvarial osteoblast:0.3,stromal cell of ovary:0.3):0.1,'
          '(hematopoietic cell:0.4,stromal cell:0.4):0):0.1,'
          '((pericyte cell:0.25,skin fibroblast:0.25):0.25,placental pericyte:0.5):0):0.05):0.1,osteoblast:0.65);')

# test random_100_250_runtime
random_100_250_newick = "((((((((P0:18,P75:18):5,(P178:18.5,P246:18.5):4.5):5,((P118:17,P33:17):7.5," \
                        "((P125:19.5,P18:19.5):1,P235:20.5):4):3.5):2," \
                        "(((P13:19,P133:19):3,(P162:17.5,P4:17.5):4.5):4,((P154:20,P93:20):3," \
                        "(P73:19.5,P85:19.5):3.5):3):4):1.5,((((P10:18.5,P28:18.5):6,(P117:19,P37:19):5.5):3.5," \
                        "(((P115:18.5,P48:18.5):4,(P138:17.5,P160:17.5):5):3,(P186:19.5,(P249:19,P47:19):0.5):6):2.5):2," \
                        "((((P101:17.5,P137:17.5):6.5,(P108:17.5,P89:17.5):6.5):2,((P155:18.5,P54:18.5):3,P165:21.5):4.5):3," \
                        "(((P119:20,P70:20):3.5,(P236:15,P61:15):8.5):3,((P187:18.5,P36:18.5):5,(P27:18,P80:18):5.5):3):2.5):1):1.5):1," \
                        "(((((P100:19.5,P157:19.5):3,(P145:18.5,P19:18.5):4):3,((P179:17.5,P42:17.5):5," \
                        "(P245:17,P52:17):5.5):3):3.5,(((P14:18.5,P144:18.5):5,((P17:20,P202:20):1,P46:21):2.5):3.5," \
                        "((P159:18,P22:18):5.5,(P244:22,P51:22):1.5):3.5):2):2,((((P111:16,P128:16):7,(P20:16,P229:16):7):1," \
                        "((P114:18,P168:18):3,P205:21):3):5.5,((((P123:14.5,P161:14.5):8.5,(P230:18,P88:18):5):4.5," \
                        "((P129:19,P32:19):3.5,(P53:16.5,P6:16.5):6):5):1,(((P140:19.5,P141:19.5):1,P234:20.5):4.5," \
                        "(P30:20,P49:20):5):3.5):1):1.5):1.5):1.5,((((((P1:19.5,P237:19.5):3,(P189:19.5,P228:19.5):3):4," \
                        "((P122:18.5,P217:18.5):3,(P57:20.5,P84:20.5):1):5):2,((((P2:17,P226:17):1.5,P69:18.5):4.5," \
                        "(P240:17.5,P241:17.5):5.5):3.5,(P81:20.5,(P90:18,P99:18):2.5):6):2):2,(((P105:17,P181:17):5.5," \
                        "(P150:17,P77:17):5.5):4.5,((P188:18,P225:18):7,(P201:18.5,P95:18.5):6.5):2):3.5):1," \
                        "(((((P11:19,P142:19):3.5,(P210:18.5,P34:18.5):4):3,(P164:16.5,P86:16.5):9):2," \
                        "(((P158:20,P180:20):2,P182:22):2,(P247:18.5,P44:18.5):5.5):3.5):2,(((P120:19.5,P16:19.5):5.5," \
                        "((P127:17.5,P64:17.5):2,P167:19.5):5.5):2.5,(((P131:18,P166:18):4,(P215:20.5,P222:20.5):1.5):2," \
                        "(P172:19,P183:19):5):3.5):2):2):2.5):1,(((((((P102:18.5,P191:18.5):3,(P132:17.5,P63:17.5):4):2.5," \
                        "(P106:19.5,(P194:16,P72:16):3.5):4.5):4,(((P112:17.5,P24:17.5):6,(P151:17.5,P71:17.5):6):0.5," \
                        "(P146:19.5,P232:19.5):4.5):4):2,((((P110:18.5,P134:18.5):5.5,((P15:17.5,P203:17.5):4.5," \
                        "(P238:17.5,P40:17.5):4.5):2):2.5,((P163:17.5,P198:17.5):5.5,(P169:19,P26:19):4):3.5):2.5," \
                        "(((P116:18,P25:18):4.5,(P62:18.5,P83:18.5):4):3,((P121:21,P45:21):3," \
                        "(P184:18.5,P39:18.5):5.5):1.5):3.5):1):2.5,(((((P103:17,P218:17):4.5,(P200:18,P56:18):3.5):5," \
                        "(((P147:18.5,P176:18.5):5,(P177:18,P66:18):5.5):2,((P209:19,P219:19):1,P223:20):5.5):1):2.5," \
                        "((((P124:18.5,P96:18.5):1,P21:19.5):5,(P185:19,P211:19):5.5):2.5,((P152:18.5,P204:18.5):4," \
                        "(P29:16,P60:16):6.5):4.5):2):1.5,((((P107:21,(P5:16.5,P67:16.5):4.5):5.5,((P199:18,P231:18):5.5," \
                        "(P92:18.5,P94:18.5):5):3):1,((P213:19,P68:19):5,(P227:20.5,(P248:18,P3:18):2.5):3.5):3.5):1.5," \
                        "(((P149:19.5,P197:19.5):2.5,(P192:19,P91:19):3):3,(P156:21,(P221:19,P74:19):2):4):4):1.5):2):1," \
                        "((((((P104:18,P78:18):5,(P113:19.5,P87:19.5):3.5):3.5,((P208:18,P239:18):5.5," \
                        "(P58:20,P98:20):3.5):3):2.5,(((P171:17,P243:17):5.5,(P216:17,P23:17):5.5):4," \
                        "((P242:18,P31:18):4,(P59:18.5,P76:18.5):3.5):4.5):2.5):1.5,((((P12:17.5,P38:17.5):4.5,P65:22):2.5," \
                        "((P135:20.5,P8:20.5):2,(P173:17,P174:17):5.5):2):4,(((P170:16,P206:16):6,(P233:18,P82:18):4):3.5," \
                        "((P207:15,P50:15):6.5,P43:21.5):4):3):2):0.5,(((((P109:18.5,P136:18.5):5,(P143:18.5,P212:18.5):5):2," \
                        "(P153:19.5,P214:19.5):6):1.5,((P126:19.5,(P130:18.5,P195:18.5):1):6,((P190:18.5,P9:18.5):5," \
                        "(P193:18,P220:18):5.5):2):1.5):2,(((P139:20.5,P224:20.5):4,((P175:15.5,P35:15.5):5.5,P55:21):3.5):3.5," \
                        "((P148:17.5,P196:17.5):7.5,((P41:18,P7:18):5.5,(P79:19,P97:19):4.5):1.5):3):1):2):2.5):1.5);"
test_statement = 'test_case("random_100_250", "complete", manhattan_distance, random_100_250_newick)'
random_100_250_runtime = timeit.timeit(test_statement, number=1, globals=globals())
assert random_100_250_runtime < 12, "your cluster_bottom_up implementation is too inefficient"
print("SUCCESS: random_100_250_runtime test passed")
