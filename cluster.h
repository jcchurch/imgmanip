#ifndef __JAMES_CLUSTER
#define __JAMES_CLUSTER
#include <stdlib.h>

struct Node {
    int isLeaf;
    double data;
    struct Node *next;
    struct Node *left;
    struct Node *right;
};

struct Cluster {
    int size;
    struct Node* head;
};

struct Node* newNode (double d) {

    struct Node *n = (struct Node*)malloc(sizeof(struct Node));

    if ( n == 0 ) {
        fprintf ( stderr, "ERROR: Failure to allocated memory for node.\n" );
        exit (1);
    }

    n->data   = d;
    n->next   = 0;
    n->left   = 0;
    n->right  = 0;
    n->isLeaf = 0;

    return n;
}

void Node_free (struct Node *n) {

    if (n != 0) {
        if ( n->left != 0 )
            Node_free( n->left );

        if ( n->right != 0 )
            Node_free( n->right );

        if ( n->next != 0 )
            Node_free( n->next );

        free(n);
        n = 0;
    }

    return;
}

int Cluster_size(struct Cluster *c) {
    return c->size;
}

void Cluster_free(struct Cluster *c) {
    if (c != 0) {
        Node_free(c->head);
        free(c);
        c = 0;
    }
    return;
}

void preorder(struct Node *n, int level) {
    int i;

    for (i = 0; i < level; i++)
        printf("    ");
    printf("%.2f\n", n->data);
    if (n->left)  preorder(n->left, level+1);
    if (n->right) preorder(n->right, level+1);
    return;
}

void traverse(struct Cluster *c) {
    if (c && c->head)
        preorder(c->head, 0);
}

int Node_count_leafs(struct Node *n) {

    if (n == 0)
        return 0;

    if (n->isLeaf)
        return 1;

    return Node_count_leafs(n->left) + Node_count_leafs(n->right);
}

double Node_sum_leafs(struct Node *n) {

    if (n == 0)
        return 0.0;

    if (n->isLeaf)
        return n->data;

    return Node_sum_leafs(n->left) + Node_sum_leafs(n->right);
}

void Cluster_merge( struct Cluster *c, int pos ) {
    int i;
    struct Node *n = newNode(0);
    c->size--;
    if (pos == 0) {
        n->right = c->head->next;
        n->left  = c->head;
        n->next  = c->head->next->next;
        c->head  = n;
    }
    else {
        struct Node *nptr = c->head;
        for (i = 0; i < (pos-1); i++)
            nptr = nptr->next;

        n->right   = nptr->next->next;
        n->left    = nptr->next;
        n->next    = nptr->next->next->next;
        nptr->next = n;
    }

    n->left->next  = 0;
    n->right->next = 0;
    n->data = Node_sum_leafs(n)/(double)Node_count_leafs(n);
    return;
}

void Cluster_put ( struct Cluster *c, double d ) {

    struct Node *n = newNode(d);
    n->isLeaf = 1;


    if ( Cluster_size(c) == 0 || d < c->head->data ) {
        n->next = c->head;
        c->head = n;
    }
    else {
        struct Node *nptr = c->head;

        while (1) {
            if (nptr->next == 0 || d < nptr->next->data)
                break;
            nptr = nptr->next;
        }

        n->next    = nptr->next;
        nptr->next = n;
    }

    c->size++;
    return;
}

struct Cluster* newCluster() {
    int i;
    struct Cluster *c = (struct Cluster*)malloc( sizeof(struct Cluster) );

    if (c == 0) {
        fprintf(stderr, "Error: failed to allocate memory for cluster.\n");
        exit(1);
    }

    c->head = 0;
    c->size = 0;
    return c;
}

void Cluster_addArray(struct Cluster *c, double *a, int n) {
    int i;
    for (i = 0; i < n; i++)
        Cluster_put(c, a[i]);
    return;
}

void Cluster_compute(struct Cluster *c) {
    int i;
    /* While more than one element in cluster list */
    while (Cluster_size(c) > 1) {
        int    small_pos   = 0;
        double small_diff;
        double diff;
        struct Node *nptr = c->head;

        /* Compute differences between each node */
        for (i = 0; i < c->size - 1; i++) {
            diff = nptr->next->data - nptr->data;
            if (i == 0) {
                small_pos  = 0;
                small_diff = diff;
            }
            else {
                if (diff < small_diff) {
                    small_pos = i;
                    small_diff = diff;
                }
            }
            nptr = nptr->next;
        }

        /* Merge two nodes with smallest difference. */
        Cluster_merge(c, small_pos);
    }
    return;
}

#endif
