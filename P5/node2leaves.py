from tree_nodes import Node, newick2nodes


def node2leaves(i_node: int, nodes: list) -> list[str]:
    """This function should call itself in a recursive way. The idea is to start
    exploring left then right. Takes as input the node id to start the recursion
    and a dictionary of nodes. Returns a list of leaf names.

    >>> node2leaves(0, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['F', 'E', 'D', 'C', 'B', 'A']
    >>> node2leaves(1, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['B', 'A']
    >>> node2leaves(2, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['A']
    >>> node2leaves(4, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['F', 'E', 'D', 'C']
    >>> node2leaves(7, newick2nodes("((A,B),((C,(D,E)),F));"))
    ['E', 'D']
    """

    # fill in the code here
    if nodes[i_node].left == 0 or nodes[i_node].right == 0:
        return [nodes[i_node].name]
    else:
        right_leaves = node2leaves(nodes[i_node].right, nodes)
        left_leaves = node2leaves(nodes[i_node].left, nodes )
    return right_leaves + left_leaves


if __name__ == "__main__":
    nodes = newick2nodes("((A,B),((C,(D,E)),F));")
    for i in range(len(nodes)):
        leaves = node2leaves(i, nodes)
        print(f"node: {i}", leaves)