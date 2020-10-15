from anytree import NodeMixin, RenderTree, Node
class Add_node(NodeMixin):  # Add Node feature
    def __init__(self, name, Combi, dEnergy, RD_sum =0, RC_sum=0,  RD_acc = 0, R_th=0, Rcsm_noD2D=0, parent=None, children=None):
        super(Add_node, self).__init__()
        self.name = name
        self.RD_sum = RD_sum
        self.RD_acc = RD_acc
        self.Combi = Combi
        self.RC_sum = RC_sum
        self.dEnergy = dEnergy
        self.parent = parent
        self.R_th = R_th
        self.Rcsm_noD2D = Rcsm_noD2D
        if children:
            self.children = children

