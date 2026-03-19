from tree_nodes import Node, newick2nodes


def leaves(nodes):
    """
    Returns a list of leaf nodes. Each leaf node is represented as a tuple: node id and node name.

    leaves(newick2nodes("((A,B),((C,(D,E)),F));"))
    [(2, 'A'), (3, 'B'), (6, 'C'), (8, 'D'), (9, 'E'), (10, 'F')]
    sisters(newick2nodes("((A,B),(((C,D),E),F));"))
    [(2, 3, 1), (7, 8, 6)]
    leaves(newick2nodes("(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"))
    [(3, 'hmgb_chite'), (4, 'hmgl_wheat'), (5, 'hmgl_trybr'), (6, 'hmgt_mouse')]
    """
    leaves = []
    for number, node in nodes.items():
        if node.name != '':
            leaves.append((number, node.name))
    return leaves


if __name__ == "__main__":
    print(leaves(newick2nodes("((A,B),((C,(D,E)),F));")))
    print(
        leaves(
            newick2nodes(
                "(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"
            )
        )
    )