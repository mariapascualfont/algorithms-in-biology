from tree_nodes import Node, newick2nodes


def sisters(nodes):
    """
    Returns a list of sister leaves. Each sister leaf is represented as a tuple: left child, right child, and parent node ids.

    sisters(newick2nodes("((A,B),((C,(D,E)),F));"))
    [(2, 3, 1), (8, 9, 7)]
    sisters(newick2nodes("((A,B),(((C,D),E),F));"))
    [(2, 3, 1), (7, 8, 6)]
    sisters(newick2nodes("(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"))
    [(3, 4, 2)]
    """

    output = []
    for number, node in nodes.items():
        left = node.left
        right = node.right
        if nodes[left].name != '' and nodes[right].name != '':
            output.append((left, right, number))

    return output


if __name__ == "__main__":
    print(sisters(newick2nodes("((A,B),((C,(D,E)),F));")))
    print(
        sisters(
            newick2nodes(
                "(((hmgb_chite:0.10,hmgl_wheat:0.25):0.20,hmgl_trybr:0.60):0.25,hmgt_mouse:0.35);"
            )
        )
    )