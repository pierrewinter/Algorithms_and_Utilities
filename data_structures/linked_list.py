# Define linked list class, then define a linked list from a list, then traverse linked list to print each value

class Node:
    def __init__(self, val):
        self.val = val
        self.next = None

def create_ll(arr):
    head = Node(arr[0])
    node = head
    for item in arr[1:]:
        node.next = Node(item)
        node = node.next
    return head




input_array = [1,5,6,3,9,8]
ll = create_ll(input_array)

head = ll
cur_node = head
while cur_node.next is not None:
    print(cur_node.val)
    cur_node = cur_node.next