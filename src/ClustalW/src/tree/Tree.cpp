/**
 * Author: Mark Larkin
 * 
 * Copyright (c) 2007 Des Higgins, Julie Thompson and Toby Gibson.  
 */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#endif
#include <sstream>
#include "Tree.h"

namespace clustalw
{

/**
 * 
 * @param firstSeq 
 * @param lastSeq 
 * @param sweight 
 */
 
void Tree::calcSeqWeights(int firstSeq, int lastSeq, vector<int>* sweight)
{
    if((int)sweight->size() < lastSeq - 1)
    {
        sweight->resize(lastSeq - 1);
    }
    
    int i, _nSeqs;
    int temp, sum;
    int *weight;
    /*
     * If there are more than three sequences....
     */
    _nSeqs = lastSeq - firstSeq;
    if ((_nSeqs >= 2) && (clustalw::userParameters->getDistanceTree() == true) && 
        (clustalw::userParameters->getNoWeights() == false))
    {
        /*
         * Calculate sequence weights based on Phylip tree.
         */

        weight = new int[lastSeq + 1];
        
        for (i = firstSeq; i < lastSeq; i++)
        {
            weight[i] = calcWeight(i);
        }

        /*
         * Normalise the weights, such that the sum of the weights = INT_SCALE_FACTOR
         */

        sum = 0;
        for (i = firstSeq; i < lastSeq; i++)
        {
            sum += weight[i];
        }

        if (sum == 0)
        {
            for (i = firstSeq; i < lastSeq; i++)
            {
                weight[i] = 1;
            }
            sum = i;
        }

        for (i = firstSeq; i < lastSeq; i++)
        {
            (*sweight)[i] = (weight[i] * INT_SCALE_FACTOR) / sum;
            if ((*sweight)[i] < 1)
            {
                (*sweight)[i] = 1;
            }
        }
        
        delete []weight;
        weight = NULL;
    }

    else
    {
        /*
         * Otherwise, use identity weights.
         */
        temp = INT_SCALE_FACTOR / _nSeqs;
        // AW 2009-05-28: goes wrong if we have more than
        // INT_SCALE_FACTOR seqs. if so, set to 1, just as above
        if (temp < 1)
            temp = 1;

        for (i = firstSeq; i < lastSeq; i++)
        {
            (*sweight)[i] = temp;
        }
    }

}

/**
 * 
 * @param alignPtr 
 * @param treeFileName 
 * @param firstSeq 
 * @param lastSeq 
 * @return 
 */
int Tree::readTree(clustalw::Alignment* alignPtr, const string& treeFileName, int firstSeq, int lastSeq)
{

    if(alignPtr == NULL || firstSeq < 0 || lastSeq < 1)
    {
        return 0;
    }
    char c;
    string name1 = "", name2 = "";
    int i, j, k;
    bool found;

    numSeq = 0;
    nnodes = 0;
    ntotal = 0;
    rootedTree = true;

    // Need to check what happens if I try to open a file that doesnt exist!
    if(treeFileName.empty())
    {
        return 0; // Cannot open, empty string!
    }
    else
    {
        file.open(treeFileName.c_str(), ifstream::in);
        if(!file.is_open())
        {
            clustalw::utilityObject->error("Cannot open output file [%s]\n", treeFileName.c_str());
        }
    }

    skipSpace(&file);
    charFromFile = file.get();
    
    if (charFromFile != '(')
    {
        clustalw::utilityObject->error("Wrong format in tree file %s\n", treeFileName.c_str());
        return 0;
    }
    file.seekg(0, ios::beg); // Put get pointer back to the begining!

    clustalw::userParameters->setDistanceTree(true);

    
    // Allocate memory for tree
    
    // Allocate memory.
    nptr = new TreeNode*[3 * (lastSeq - firstSeq + 1)];
    ptrs = new TreeNode*[3 * (lastSeq - firstSeq + 1)];
    lptr = new TreeNode*[lastSeq - firstSeq + 1];
    olptr = new TreeNode*[lastSeq + 1];
    
    seqTree = avail();
    setInfo(seqTree, NULL, 0, string(""), 0.0);

    createTree(seqTree, NULL, &file);
    file.close();

    if (numSeq != lastSeq - firstSeq)
    {
        stringstream ss;
        ss << "tree not compatible with alignment\n(" << lastSeq - firstSeq 
        << " sequences in alignment and "<< numSeq <<" in tree\n";
        string errorMes;
        ss >> errorMes;
        clustalw::utilityObject->error(errorMes.c_str());
        return 0;
    }
   
    // If the tree is unrooted, reroot the tree - ie. minimise the difference
    // between the mean root->leaf distances for the left and right branches of
    // the tree.     

    if (clustalw::userParameters->getDistanceTree() == false)
    {
        if (rootedTree == false)
        {
            clustalw::utilityObject->error("input tree is unrooted and has no distances.\nCannot align sequences");
            return 0;
        }
    }

    if (rootedTree == false)
    {
        root = reRoot(seqTree, lastSeq - firstSeq + 1);
    }
    else
    {
        root = seqTree;
    }
   
    // calculate the 'order' of each node.
    int nameLength;
    string nameSeq;
    orderNodes();

    if (numSeq >= 2)
    {
        // If there are more than three sequences....
        // assign the sequence nodes (in the same order as in the alignment file)
        for (i = firstSeq; i < lastSeq; i++)
        {
            nameLength = alignPtr->getName(i + 1).length();
            nameSeq = alignPtr->getName(i + 1);
            name1 = "";
            name2 = "";
            if (nameLength > clustalw::MAXNAMES)
            {
                stringstream ss;
                ss << "name " << nameSeq << " is too long for PHYLIP tree format (max " 
                   << clustalw::MAXNAMES << " chars)";
                string msg; 
                ss >> msg; 
                clustalw::utilityObject->warning(msg.c_str());// Mark change 17-5-07
            }

            for (k = 0; k < nameLength && k < clustalw::MAXNAMES; k++)
            {
                c = nameSeq[k];
                if ((c > 0x40) && (c < 0x5b)) 
                {
                    c = c | 0x20;
                }
                if (c == ' ')
                {
                    c = '_';
                }
                name2 += c;
            }
            found = false;

            for (j = 0; j < numSeq; j++)
            {
                name1 = "";
                for (k = 0; k < (int)lptr[j]->name.length() && k < clustalw::MAXNAMES; k++)
                {
                    c = lptr[j]->name[k];
                    if ((c > 0x40) && (c < 0x5b))
                    {
                        c = c | 0x20;
                    }

                    name1 += c;
                }

                if (name1.compare(name2) == 0)
                {
                    olptr[i] = lptr[j];
                    found = true;
                }
            }

            if (found == false)
            {
                utilityObject->error("tree not compatible with alignment:\n %s not found\n", name2.c_str());
                return 0;
            }
        }

    }

    return (1);


}

/**
 * 
 * @param firstSeq 
 * @param lastSeq 
 * @return 
 */
auto_ptr<AlignmentSteps> Tree::createSets(int firstSeq, int lastSeq)
{
    auto_ptr<AlignmentSteps> progAlignSteps;
    progAlignSteps.reset(new AlignmentSteps);
    
    int i, j, _nSeqs;

    numSets = 0;
    _nSeqs = lastSeq - firstSeq;
    if (_nSeqs >= 2)
    {
        // If there are more than three sequences....
        groups = new int[_nSeqs + 1];
        groupSeqs(root, groups, _nSeqs, progAlignSteps.get());

        delete []groups;
        groups  = NULL;

    }

    else
    {
        groups = new int[_nSeqs + 1];
        for (i = 0; i < _nSeqs - 1; i++)
        {
            for (j = 0; j < _nSeqs; j++)
                if (j <= i)
                {
                    groups[j] = 1;
                }
                else if (j == i + 1)
                {
                    groups[j] = 2;
                }
                else
                {
                    groups[j] = 0;
                }

            progAlignSteps->saveSet(_nSeqs, groups);
        }
        delete []groups;
        groups  = NULL;
    }    
    
    return progAlignSteps;
}

/**
 * calcSimilarities changes the distMat.
 * @param alignPtr 
 * @param distMat 
 * @return 
 */
int Tree::calcSimilarities(clustalw::Alignment* alignPtr, clustalw::DistMatrix* distMat)
{
    int depth = 0, i, j, k, n;
    bool found;
    int nerrs, seq1[MAXERRS], seq2[MAXERRS];
    TreeNode *p,  **pathToRoot;
    float dist;
    float *distToNode, badDist[MAXERRS];
    double **dmat;
    ostringstream err1;
    char reply;
    int nSeqs = alignPtr->getNumSeqs();
    
    pathToRoot = new TreeNode*[nSeqs];
    distToNode = new float[nSeqs];
    dmat = new double*[nSeqs];

    for (i = 0; i < nSeqs; i++)
    {
        dmat[i] = new double[nSeqs];
    }

    if (nSeqs >= 2)
    {
        /*
         * for each leaf, determine all nodes between the leaf and the root;
         */
        for (i = 0; i < nSeqs; i++)
        {
            depth = 0;dist = 0.0;
            p = olptr[i];
            while (p != NULL)
            {
                pathToRoot[depth] = p;
                dist += p->dist;
                distToNode[depth] = dist;
                p = p->parent;
                depth++;
            }

            /*
             * for each pair....
             */
            for (j = 0; j < i; j++)
            {
                p = olptr[j];
                dist = 0.0;
                /*
                 * find the common ancestor.
                 */
                found = false;
                n = 0;
                while ((found == false) && (p->parent != NULL))
                {
                    for (k = 0; k < depth; k++)
                    if (p->parent == pathToRoot[k])
                    {
                        found = true;
                        n = k;
                    }

                    dist += p->dist;
                    p = p->parent;
                }

                dmat[i][j] = dist + distToNode[n - 1];
            }
        }

        nerrs = 0;
        for (i = 0; i < nSeqs; i++)
        {
            dmat[i][i] = 0.0;
            for (j = 0; j < i; j++)
            {
                if (dmat[i][j] < 0.01)
                {
                    dmat[i][j] = 0.01;
                }
                if (dmat[i][j] > 1.0)
                {
                    if (dmat[i][j] > 1.1 && nerrs < MAXERRS)
                    {
                        seq1[nerrs] = i;
                        seq2[nerrs] = j;
                        badDist[nerrs] = dmat[i][j];
                        nerrs++;
                    }
                    dmat[i][j] = 1.0;
                }
            }
        }
        if (nerrs > 0 && !userParameters->getGui())
        {
            string errMess = "The following sequences are too divergent to be aligned:\n";
            
            for (i = 0; i < nerrs && i < 5; i++)
            {
                err1 << "           " << alignPtr->getName(seq1[i] + 1) << " and " 
                     << alignPtr->getName(seq2[i] + 1) << " (distance " 
                     << setprecision(3) << badDist[i] << ")\n"; 
            }
            errMess += err1.str();
            errMess += "(All distances should be between 0.0 and 1.0)\n";
            errMess += "This may not be fatal but you have been warned!\n";
            errMess += "SUGGESTION: Remove one or more problem sequences and try again";
            
            if (clustalw::userParameters->getInteractive())
            {
                reply = clustalw::utilityObject->promptForYesNo(errMess.c_str(), "Continue ");
            }
            else
            {
                reply = 'y';
            }
            if ((reply != 'y') && (reply != 'Y'))
            {
                return 0;
            }
        }
    }
    else
    {
        for (i = 0; i < nSeqs; i++)
        {
            for (j = 0; j < i; j++)
            {
                dmat[i][j] = (*distMat)(i + 1, j + 1);
            }
        }
    }
    delete []pathToRoot;
    delete []distToNode;
    
    double value;
    for (i = 0; i < nSeqs; i++)
    {
        distMat->SetAt(i + 1, i + 1, 0.0);
        for (j = 0; j < i; j++)
        {
            value = 100.0 - (dmat[i][j]) * 100.0;
            distMat->SetAt(i + 1, j + 1, value);
            distMat->SetAt(j + 1, i + 1, value);
        }
    }

    for (i = 0; i < nSeqs; i++)
    {
        delete [] dmat[i];
    }

    delete []dmat;

    return 1;
    
}

/** *************************************************************************
 * Private functions!!!!!!!!!!!!!!!                                         *
 ****************************************************************************/
 
/**
 * 
 * @param ptree 
 * @param parent 
 * @param file 
 */
void Tree::createTree(clustalw::TreeNode* ptree, clustalw::TreeNode* parent, ifstream* file)
{
    TreeNode* p;

    int i, type;
    float dist;
    string name;

    
    // is this a node or a leaf ?
    skipSpace(file);
    charFromFile = file->get();
    if (charFromFile == '(')
    {
        // this must be a node....
        type = NODE;
        name = "";
        ptrs[ntotal] = nptr[nnodes] = ptree;
        nnodes++;
        ntotal++;

        createNode(ptree, parent);

        p = ptree->left;
        createTree(p, ptree, file);

        if (charFromFile == ',')
        {
            p = ptree->right;
            createTree(p, ptree, file);
            if (charFromFile == ',')
            {
                ptree = insertNode(ptree);
                ptrs[ntotal] = nptr[nnodes] = ptree;
                nnodes++;
                ntotal++;
                p = ptree->right;
                createTree(p, ptree, file);
                rootedTree = false;
            }
        }

        skipSpace(file);
        charFromFile = file->get();
    }
    
    // ...otherwise, this is a leaf
    
    else
    {
        type = LEAF;
        ptrs[ntotal++] = lptr[numSeq++] = ptree;        
        // get the sequence name
        name = "";
        name += charFromFile;
        charFromFile = file->get();
        
        i = 1;
        while ((charFromFile != ':') && (charFromFile != ',') && (charFromFile != ')'))
        {
            if (i < MAXNAMES)
            {
                name += charFromFile;
                i++;
            }
            charFromFile = file->get();
        }

        if (charFromFile != ':')
        {
            clustalw::userParameters->setDistanceTree(false);
            dist = 0.0;
        }
    }
    
    // get the distance information
    
    dist = 0.0;
    if (charFromFile == ':')
    {
        skipSpace(file);
        (*file) >> dist;
        skipSpace(file);
        charFromFile = file->get();
    }
    setInfo(ptree, parent, type, name, dist);
}

/**
 * 
 * @param pptr 
 * @param parent 
 */
void Tree::createNode(TreeNode* pptr, TreeNode* parent)
{
    TreeNode* t;
    pptr->parent = parent;
    t = avail();
    pptr->left = t;
    t = avail();
    pptr->right = t;
}

/**
 * 
 * @param pptr 
 * @return 
 */
TreeNode* Tree::insertNode(TreeNode* pptr)
{

    TreeNode* newnode;
    newnode = avail();
    //createNode(newnode, pptr->parent);
    //copy/paste createNode with changes - BEGIN //valgrind - duplicate/unused TreeNode?
	TreeNode* t;
	newnode->parent = pptr->parent;
	TreeNode* father;
	if (pptr->parent != NULL)
	{

		father = pptr->parent;
		if (father->right == pptr)
		{
			father->right = newnode;
		}
		else
		{
			father->left = newnode;
		}
	}
	t = avail();
	t->parent = newnode; //inserted ...
	newnode->right = t;
	//END

    newnode->left = pptr;
    pptr->parent = newnode;

    if (NODE)
    {
    	newnode->left = NULL;
    	newnode->right = NULL;
    }

    //setInfo(newnode, pptr->parent, NODE, "", 0.0);
    return newnode;
}

/**
 * 
 * @param p 
 */
void Tree::clearTree(clustalw::TreeNode* p)
{
    clearTreeNodes(p);
    delete [] nptr;
    nptr  = NULL;
    delete [] ptrs;
    ptrs  = NULL;
    delete [] lptr;
    lptr  = NULL;
    delete [] olptr;
    olptr  = NULL;
}

/**
 * 
 * @param p 
 */
void Tree::clearTreeNodes(clustalw::TreeNode* p)
{
    if (p == NULL)
    {
        p = root;
    }
    if (p->left != NULL)
    {
        clearTreeNodes(p->left);
    }
    if (p->right != NULL)
    {
        clearTreeNodes(p->right);
    }
    p->left = NULL;
    p->right = NULL;
    
    delete p;
    p  = NULL;
}




/**
 * 
 * @param 
 * @param 
 * @return 
 */
void Tree::debugPrintAllNodes(int nseqs)
{
  clustalw::TreeNode *p;
  int i;
   float diff, maxDist;
 
    cerr << "\nDEBUG: reportAllNodes\n";
    for (i = 0; i < ntotal; i++) {
        p = ptrs[i];
        //        ios::sync_with_stdio();

        // same design as TreeNode
        if (p->parent == NULL)
            diff = calcRootMean(p, &maxDist);
        else
            diff = calcMean(p, &maxDist, nseqs);
        fprintf(stdout,
                "i=%d p=%p: parent=%p left=%p right=%p dist=%f diff=%f\n",
                i, (void*)p, (void*)p->parent, (void*)p->left, (void*)p->right,
		p->dist, diff);
    }
}





/**
 * 
 * @param ptree 
 * @param nseqs 
 * @return 
 */
clustalw::TreeNode* Tree::reRoot(clustalw::TreeNode* ptree, int nseqs)
{
    clustalw::TreeNode *p, *rootNode, *rootPtr;
    float diff, minDiff = 0.0, minDepth = 1.0, maxDist;
    int i;
    bool first = true;

    // find the difference between the means of leaf->node
    // distances on the left and on the right of each node
    rootPtr = ptree;
    for (i = 0; i < ntotal; i++)
    {
        p = ptrs[i];
        if (p->parent == NULL)
        {
            /* AW Bug 94: p->parent must be chosen as rootNode
               (non-optimized executables (-O0) never do), otherwise
               insertRoot fails.
               Is p->parent == NULL valid at all?
               Simplest thing for now is to continue here. Tree code
               needs serious dismantling anyway. See debugPrintAllNodes
            */

            continue;
            //diff = calcRootMean(p, &maxDist);
        }
        else
        {
            diff = calcMean(p, &maxDist, nseqs);
        }

        if ((diff == 0) || ((diff > 0) && (diff < 2 *p->dist)))
        {
            if ((maxDist < minDepth) || (first == true))
            {
                first = false;
                rootPtr = p;
                minDepth = maxDist;
                minDiff = diff;
            }
        }

    }

    
    // insert a new node as the ancestor of the node which produces the shallowest
    // tree.
    /* AW Bug 94: could also be prevented here */
    if (rootPtr == ptree)
    {
        minDiff = rootPtr->left->dist + rootPtr->right->dist;
        rootPtr = rootPtr->right;
    }
    rootNode = insertRoot(rootPtr, minDiff);

    diff = calcRootMean(rootNode, &maxDist);

    return rootNode;
}

/**
 * 
 * @param p 
 * @param diff 
 * @return 
 */
clustalw::TreeNode* Tree::insertRoot(clustalw::TreeNode* p, float diff)
{
    clustalw::TreeNode *newp, *prev, *q, *t;
    float dist, prevDist, td;

    newp = avail();

    if (p->parent==NULL) {
        // AW bug 94: question remains if access here should be handled differently
        cerr << "\n\n*** INTERNAL ERROR: Tree::insertRoot: TreeNode p->parent is NULL\n";
        cerr << "To help us fix this bug, please send sequence file and used options to clustalw@ucd.ie\n";
        throw 1;
    }

    t = p->parent;
    prevDist = t->dist;

    p->parent = newp;

    dist = p->dist;

    p->dist = diff / 2;

    if (p->dist < 0.0)
    {
        p->dist = 0.0;
    }

    if (p->dist > dist)
    {
        p->dist = dist;
    }

    t->dist = dist - p->dist;

    newp->left = t;
    newp->right = p;
    newp->parent = NULL;
    newp->dist = 0.0;
    newp->leaf = NODE;

    if (t->left == p)
    {
        t->left = t->parent;
    }
    else
    {
        t->right = t->parent;
    }

    prev = t;
    q = t->parent;

    t->parent = newp;

    while (q != NULL)
    {
        if (q->left == prev)
        {
            q->left = q->parent;
            q->parent = prev;
            td = q->dist;
            q->dist = prevDist;
            prevDist = td;
            prev = q;
            q = q->left;
        }
        else
        {
            q->right = q->parent;
            q->parent = prev;
            td = q->dist;
            q->dist = prevDist;
            prevDist = td;
            prev = q;
            q = q->right;
        }
    }

    /*
     * remove the old root node
     */
    q = prev;
    if (q->left == NULL)
    {
        dist = q->dist;
        q = q->right;
        q->dist += dist;
        q->parent = prev->parent;

        if (prev->parent->left == prev)
        {
            prev->parent->left = q;
        }
        else
        {
            prev->parent->right = q;
        }

        prev->right = NULL;
    }
    else
    {
        dist = q->dist;
        q = q->left;
        q->dist += dist;
        q->parent = prev->parent;

        if (prev->parent->left == prev)
        {
            prev->parent->left = q;
        }
        else
        {
            prev->parent->right = q;
        }

        prev->left = NULL;
    }

    return (newp);
}

/**
 * 
 * @param root 
 * @param maxDist 
 * @return 
 */
float Tree::calcRootMean(clustalw::TreeNode* root, float *maxDist)
{
    float dist, leftSum = 0.0, rightSum = 0.0, leftMean, rightMean, diff;
    clustalw::TreeNode* p;
    int i;
    int numLeft, numRight;
    int direction;
    
    // for each leaf, determine whether the leaf is left or right of the root.

    dist = (*maxDist) = 0;
    numLeft = numRight = 0;
    for (i = 0; i < numSeq; i++)
    {
        p = lptr[i];
        dist = 0.0;
        while (p->parent != root)
        {
            dist += p->dist;
            p = p->parent;
        }

        if (p == root->left)
        {
            direction = LEFT;
        }
        else
        {
            direction = RIGHT;
        }

        dist += p->dist;

        if (direction == LEFT)
        {
            leftSum += dist;
            numLeft++;
        }
        else
        {
            rightSum += dist;
            numRight++;
        }
        if (dist > (*maxDist))
        {
            *maxDist = dist;
        }
    }

    leftMean = leftSum / numLeft;
    rightMean = rightSum / numRight;

    diff = leftMean - rightMean;
    return diff;
}

/**
 * 
 * @param nptr 
 * @param maxDist 
 * @param nSeqs 
 * @return 
 */
float Tree::calcMean(clustalw::TreeNode* nptr, float *maxDist, int nSeqs)
{
    float dist, leftSum = 0.0, rightSum = 0.0, leftMean, rightMean, diff;
    clustalw::TreeNode* p;  
    clustalw::TreeNode** pathToRoot;
    float *distToNode;
    int depth = 0, i, j, n = 0;
    int numLeft, numRight;
    int direction;
    bool found;

    pathToRoot = new clustalw::TreeNode*[nSeqs];
    distToNode = new float[nSeqs];
    // determine all nodes between the selected node and the root;
  
    depth = 0;
    (*maxDist) = dist = 0.0;
    numLeft = numRight = 0;
    p = nptr;
    while (p != NULL)
    {
        pathToRoot[depth] = p;
        dist += p->dist;
        distToNode[depth] = dist;
        p = p->parent;
        depth++;
    }

    // for each leaf, determine whether the leaf is left or right of the node.
    // (RIGHT = descendant, LEFT = not descendant)
 
    for (i = 0; i < numSeq; i++)
    {
        p = lptr[i];
        if (p == nptr)
        {
            direction = RIGHT;
            dist = 0.0;
        }
        else
        {
            direction = LEFT;
            dist = 0.0;
            
            // find the common ancestor. 
            
            found = false;
            n = 0;
            while ((found == false) && (p->parent != NULL))
            {
                for (j = 0; j < depth; j++)
                if (p->parent == pathToRoot[j])
                {
                    found = true;
                    n = j;
                }

                dist += p->dist;
                p = p->parent;
            }
            if (p == nptr)
            {
                direction = RIGHT;
            }
        }

        if (direction == LEFT)
        {
            leftSum += dist;
            leftSum += distToNode[n - 1];
            numLeft++;
        }
        else
        {
            rightSum += dist;
            numRight++;
        }

        if (dist > (*maxDist))
        {
            *maxDist = dist;
        }
    }

    delete [] distToNode;
    distToNode = NULL;
    delete [] pathToRoot;
    pathToRoot = NULL;
    
    leftMean = leftSum / numLeft;
    rightMean = rightSum / numRight;

    diff = leftMean - rightMean;
    return (diff);
}

/**
 * 
 */
void Tree::orderNodes()
{
    int i;
    clustalw::TreeNode* p;

    for (i = 0; i < numSeq; i++)
    {
        p = lptr[i];
        while (p != NULL)
        {
            p->order++;
            p = p->parent;
        }
    }
}

/**
 * 
 * @param leaf 
 * @return 
 */
int Tree::calcWeight(int leaf)
{
    clustalw::TreeNode* p;
    float weight = 0.0;

    p = olptr[leaf];
    while (p->parent != NULL)
    {
        weight += p->dist / p->order;
        p = p->parent;
    }

    weight *= 100.0;

    return ((int)weight);
}

/**
 * skipSpace is used to skip all the spaces at the begining of a file. The next read
 * will give a character other than a space.
 * @param file 
 */
void Tree::skipSpace(ifstream* file)
{
    char c;

    do
    {
        c = file->get();
    }
    while (isspace(c));

    file->putback(c);
}

/**
 * 
 * @param p 
 * @param nextGroups 
 * @param nSeqs 
 * @param stepsPtr 
 */
void Tree::groupSeqs(clustalw::TreeNode* p, int *nextGroups, int nSeqs, AlignmentSteps* stepsPtr)
{
    int i;
    int *tmpGroups;

    tmpGroups = new int[nSeqs + 1];
    for (i = 0; i < nSeqs; i++)
    {
        tmpGroups[i] = 0;
    }

    if (p->left != NULL)
    {
        if (p->left->leaf == NODE)
        {
            groupSeqs(p->left, nextGroups, nSeqs, stepsPtr);

            for (i = 0; i < nSeqs; i++)
                if (nextGroups[i] != 0)
                {
                    tmpGroups[i] = 1;
                }
        }
        else
        {
            markGroup1(p->left, tmpGroups, nSeqs);
        }

    }

    if (p->right != NULL)
    {
        if (p->right->leaf == NODE)
        {
            groupSeqs(p->right, nextGroups, nSeqs, stepsPtr);
            for (i = 0; i < nSeqs; i++)
                if (nextGroups[i] != 0)
                {
                    tmpGroups[i] = 2;
                }
        }
        else
        {
            markGroup2(p->right, tmpGroups, nSeqs);
        }
        stepsPtr->saveSet(nSeqs, tmpGroups);
    }

    for (i = 0; i < nSeqs; i++)
    {
        nextGroups[i] = tmpGroups[i];
    }

    delete [] tmpGroups;
    tmpGroups  = NULL;
}

/**
 * 
 * @param p 
 * @param groups 
 * @param n 
 */
void Tree::markGroup1(clustalw::TreeNode* p, int *groups, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (olptr[i] == p)
        {
            groups[i] = 1;
        }
        else
        {
            groups[i] = 0;
        }
    }
}

/**
 * 
 * @param p 
 * @param groups 
 * @param n 
 */
void Tree::markGroup2(clustalw::TreeNode* p, int *groups, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (olptr[i] == p)
        {
            groups[i] = 2;
        }
        else if (groups[i] != 0)
        {
            groups[i] = 1;
        }
    }
}

/**
 * 
 * @return 
 */
clustalw::TreeNode* Tree::avail()
{
    clustalw::TreeNode* p;
    p = new TreeNode; 
    p->left = NULL;
    p->right = NULL;
    p->parent = NULL;
    p->dist = 0.0;
    p->leaf = 0;
    p->order = 0;
    p->name = "";
    return (p);
}

/**
 * 
 * @param p 
 * @param parent 
 * @param pleaf 
 * @param pname 
 * @param pdist 
 */
void Tree::setInfo(TreeNode* p, TreeNode* parent, int pleaf, string pname, float
                     pdist)
{
    p->parent = parent;
    p->leaf = pleaf;
    p->dist = pdist;
    p->order = 0;
    p->name = pname;
    if (p->leaf == true)
    {
        p->left = NULL;
        p->right = NULL;
    }
}
          
}
