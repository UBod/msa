#include "muscle.h"
#include "tree.h"
#include <math.h>

#define TRACE 0

/***
Node has 0 to 3 neighbors:
	0 neighbors:	singleton root
	1 neighbor:		leaf, neighbor is parent
	2 neigbors:		non-singleton root
	3 neighbors:	internal node (other than root)

Minimal rooted tree is single node.
Minimal unrooted tree is single edge.
Leaf node always has nulls in neighbors 2 and 3, neighbor 1 is parent.
When tree is rooted, neighbor 1=parent, 2=left, 3=right.
***/

void Tree::AssertAreNeighbors(unsigned uNodeIndex1, unsigned uNodeIndex2) const
	{
	if (uNodeIndex1 >= m_uNodeCount || uNodeIndex2 >= m_uNodeCount)
		Quit("AssertAreNeighbors(%u,%u), are %u nodes",
		  uNodeIndex1, uNodeIndex2, m_uNodeCount);

	if (m_uNeighbor1[uNodeIndex1] != uNodeIndex2 &&
	  m_uNeighbor2[uNodeIndex1] != uNodeIndex2 &&
	  m_uNeighbor3[uNodeIndex1] != uNodeIndex2)
		{
		LogMe();
		Quit("AssertAreNeighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
		}

	if (m_uNeighbor1[uNodeIndex2] != uNodeIndex1 &&
	  m_uNeighbor2[uNodeIndex2] != uNodeIndex1 &&
	  m_uNeighbor3[uNodeIndex2] != uNodeIndex1)
		{
		LogMe();
		Quit("AssertAreNeighbors(%u,%u) failed", uNodeIndex1, uNodeIndex2);
		}

	bool Has12 = HasEdgeLength(uNodeIndex1, uNodeIndex2);
	bool Has21 = HasEdgeLength(uNodeIndex2, uNodeIndex1);
	if (Has12 != Has21)
		{
		HasEdgeLength(uNodeIndex1, uNodeIndex2);
		HasEdgeLength(uNodeIndex2, uNodeIndex1);
		LogMe();
		Log("HasEdgeLength(%u, %u)=%c HasEdgeLength(%u, %u)=%c\n",
		  uNodeIndex1,
		  uNodeIndex2,
		  Has12 ? 'T' : 'F',
		  uNodeIndex2,
		  uNodeIndex1,
		  Has21 ? 'T' : 'F');

		Quit("Tree::AssertAreNeighbors, HasEdgeLength not symmetric");
		}

	if (Has12)
		{
		double d12 = GetEdgeLength(uNodeIndex1, uNodeIndex2);
		double d21 = GetEdgeLength(uNodeIndex2, uNodeIndex1);
		if (d12 != d21)
			{
			LogMe();
			Quit("Tree::AssertAreNeighbors, Edge length disagrees %u-%u=%.3g, %u-%u=%.3g",
			  uNodeIndex1, uNodeIndex2, d12,
			  uNodeIndex2, uNodeIndex1, d21);
			}
		}
	}

void Tree::ValidateNode(unsigned uNodeIndex) const
	{
	if (uNodeIndex >= m_uNodeCount)
		Quit("ValidateNode(%u), %u nodes", uNodeIndex, m_uNodeCount);

	const unsigned uNeighborCount = GetNeighborCount(uNodeIndex);

	if (2 == uNeighborCount)
		{
		if (!m_bRooted)
			{
			LogMe();
			Quit("Tree::ValidateNode: Node %u has two neighbors, tree is not rooted",
			 uNodeIndex);
			}
		if (uNodeIndex != m_uRootNodeIndex)
			{
			LogMe();
			Quit("Tree::ValidateNode: Node %u has two neighbors, but not root node=%u",
			 uNodeIndex, m_uRootNodeIndex);
			}
		}

	const unsigned n1 = m_uNeighbor1[uNodeIndex];
	const unsigned n2 = m_uNeighbor2[uNodeIndex];
	const unsigned n3 = m_uNeighbor3[uNodeIndex];

	if (NULL_NEIGHBOR == n2 && NULL_NEIGHBOR != n3)
		{
		LogMe();
		Quit("Tree::ValidateNode, n2=null, n3!=null", uNodeIndex);
		}
	if (NULL_NEIGHBOR == n3 && NULL_NEIGHBOR != n2)
		{
		LogMe();
		Quit("Tree::ValidateNode, n3=null, n2!=null", uNodeIndex);
		}

	if (n1 != NULL_NEIGHBOR)
		AssertAreNeighbors(uNodeIndex, n1);
	if (n2 != NULL_NEIGHBOR)
		AssertAreNeighbors(uNodeIndex, n2);
	if (n3 != NULL_NEIGHBOR)
		AssertAreNeighbors(uNodeIndex, n3);

	if (n1 != NULL_NEIGHBOR && (n1 == n2 || n1 == n3))
		{
		LogMe();
		Quit("Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
		}
	if (n2 != NULL_NEIGHBOR && (n2 == n1 || n2 == n3))
		{
		LogMe();
		Quit("Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
		}
	if (n3 != NULL_NEIGHBOR && (n3 == n1 || n3 == n2))
		{
		LogMe();
		Quit("Tree::ValidateNode, duplicate neighbors in node %u", uNodeIndex);
		}

	if (IsRooted())
		{
		if (NULL_NEIGHBOR == GetParent(uNodeIndex))
			{
			if (uNodeIndex != m_uRootNodeIndex)
				{
				LogMe();
				Quit("Tree::ValiateNode(%u), no parent", uNodeIndex);
				}
			}
		else if (GetLeft(GetParent(uNodeIndex)) != uNodeIndex &&
		  GetRight(GetParent(uNodeIndex)) != uNodeIndex)
			{
			LogMe();
			Quit("Tree::ValidateNode(%u), parent / child mismatch", uNodeIndex);
			}
		}
	}

void Tree::Validate() const
	{
	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		ValidateNode(uNodeIndex);
	}

bool Tree::IsEdge(unsigned uNodeIndex1, unsigned uNodeIndex2) const
	{
	assert(uNodeIndex1 < m_uNodeCount && uNodeIndex2 < m_uNodeCount);

	return m_uNeighbor1[uNodeIndex1] == uNodeIndex2 ||
	  m_uNeighbor2[uNodeIndex1] == uNodeIndex2 ||
	  m_uNeighbor3[uNodeIndex1] == uNodeIndex2;
	}

double Tree::GetEdgeLength(unsigned uNodeIndex1, unsigned uNodeIndex2) const
	{
	assert(uNodeIndex1 < m_uNodeCount && uNodeIndex2 < m_uNodeCount);
	if (!HasEdgeLength(uNodeIndex1, uNodeIndex2))
		{
		LogMe();
		Quit("Missing edge length in tree %u-%u", uNodeIndex1, uNodeIndex2);
		}

	if (m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
		return m_dEdgeLength1[uNodeIndex1];
	else if (m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
		return m_dEdgeLength2[uNodeIndex1];
	assert(m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
	return m_dEdgeLength3[uNodeIndex1];
	}

void Tree::ExpandCache()
	{
	const unsigned uNodeCount = 100;
	unsigned uNewCacheCount = m_uCacheCount + uNodeCount;
	unsigned *uNewNeighbor1 = new unsigned[uNewCacheCount];
	unsigned *uNewNeighbor2 = new unsigned[uNewCacheCount];
	unsigned *uNewNeighbor3 = new unsigned[uNewCacheCount];

	unsigned *uNewIds = new unsigned[uNewCacheCount];
	memset(uNewIds, 0xff, uNewCacheCount*sizeof(unsigned));

	double *dNewEdgeLength1 = new double[uNewCacheCount];
	double *dNewEdgeLength2 = new double[uNewCacheCount];
	double *dNewEdgeLength3 = new double[uNewCacheCount];
	double *dNewHeight = new double[uNewCacheCount];

	bool *bNewHasEdgeLength1 = new bool[uNewCacheCount];
	bool *bNewHasEdgeLength2 = new bool[uNewCacheCount];
	bool *bNewHasEdgeLength3 = new bool[uNewCacheCount];
	bool *bNewHasHeight = new bool[uNewCacheCount];

	char **ptrNewName = new char *[uNewCacheCount];
	memset(ptrNewName, 0, uNewCacheCount*sizeof(char *));

	if (m_uCacheCount > 0)
		{
		const unsigned uUnsignedBytes = m_uCacheCount*sizeof(unsigned);
		memcpy(uNewNeighbor1, m_uNeighbor1, uUnsignedBytes);
		memcpy(uNewNeighbor2, m_uNeighbor2, uUnsignedBytes);
		memcpy(uNewNeighbor3, m_uNeighbor3, uUnsignedBytes);

		memcpy(uNewIds, m_Ids, uUnsignedBytes);

		const unsigned uEdgeBytes = m_uCacheCount*sizeof(double);
		memcpy(dNewEdgeLength1, m_dEdgeLength1, uEdgeBytes);
		memcpy(dNewEdgeLength2, m_dEdgeLength2, uEdgeBytes);
		memcpy(dNewEdgeLength3, m_dEdgeLength3, uEdgeBytes);
		memcpy(dNewHeight, m_dHeight, uEdgeBytes);

		const unsigned uBoolBytes = m_uCacheCount*sizeof(bool);
		memcpy(bNewHasEdgeLength1, m_bHasEdgeLength1, uBoolBytes);
		memcpy(bNewHasEdgeLength2, m_bHasEdgeLength2, uBoolBytes);
		memcpy(bNewHasEdgeLength3, m_bHasEdgeLength3, uBoolBytes);
		memcpy(bNewHasHeight, m_bHasHeight, uBoolBytes);

		const unsigned uNameBytes = m_uCacheCount*sizeof(char *);
		memcpy(ptrNewName, m_ptrName, uNameBytes);

		delete[] m_uNeighbor1;
		delete[] m_uNeighbor2;
		delete[] m_uNeighbor3;

		delete[] m_Ids;

		delete[] m_dEdgeLength1;
		delete[] m_dEdgeLength2;
		delete[] m_dEdgeLength3;

		delete[] m_bHasEdgeLength1;
		delete[] m_bHasEdgeLength2;
		delete[] m_bHasEdgeLength3;
		delete[] m_bHasHeight;

		delete[] m_ptrName;
		}
	m_uCacheCount = uNewCacheCount;
	m_uNeighbor1 = uNewNeighbor1;
	m_uNeighbor2 = uNewNeighbor2;
	m_uNeighbor3 = uNewNeighbor3;
	m_Ids = uNewIds;
	m_dEdgeLength1 = dNewEdgeLength1;
	m_dEdgeLength2 = dNewEdgeLength2;
	m_dEdgeLength3 = dNewEdgeLength3;
	m_dHeight = dNewHeight;
	m_bHasEdgeLength1 = bNewHasEdgeLength1;
	m_bHasEdgeLength2 = bNewHasEdgeLength2;
	m_bHasEdgeLength3 = bNewHasEdgeLength3;
	m_bHasHeight = bNewHasHeight;
	m_ptrName = ptrNewName;
	}

// Creates tree with single node, no edges.
// Root node always has index 0.
void Tree::CreateRooted()
	{
	Clear();
	ExpandCache();
	m_uNodeCount = 1;

	m_uNeighbor1[0] = NULL_NEIGHBOR;
	m_uNeighbor2[0] = NULL_NEIGHBOR;
	m_uNeighbor3[0] = NULL_NEIGHBOR;

	m_bHasEdgeLength1[0] = false;
	m_bHasEdgeLength2[0] = false;
	m_bHasEdgeLength3[0] = false;
	m_bHasHeight[0] = false;

	m_uRootNodeIndex = 0;
	m_bRooted = true;

#if	DEBUG
	Validate();
#endif
	}

// Creates unrooted tree with single edge.
// Nodes for that edge are always 0 and 1.
void Tree::CreateUnrooted(double dEdgeLength)
	{
	Clear();
	ExpandCache();

	m_uNeighbor1[0] = 1;
	m_uNeighbor2[0] = NULL_NEIGHBOR;
	m_uNeighbor3[0] = NULL_NEIGHBOR;

	m_uNeighbor1[1] = 0;
	m_uNeighbor2[1] = NULL_NEIGHBOR;
	m_uNeighbor3[1] = NULL_NEIGHBOR;

	m_dEdgeLength1[0] = dEdgeLength;
	m_dEdgeLength1[1] = dEdgeLength;

	m_bHasEdgeLength1[0] = true;
	m_bHasEdgeLength1[1] = true;

	m_bRooted = false;

#if	DEBUG
	Validate();
#endif
	}

void Tree::SetLeafName(unsigned uNodeIndex, const char *ptrName)
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	free(m_ptrName[uNodeIndex]);
	m_ptrName[uNodeIndex] = strsave(ptrName);
	}

void Tree::SetLeafId(unsigned uNodeIndex, unsigned uId)
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	m_Ids[uNodeIndex] = uId;
	}

const char *Tree::GetLeafName(unsigned uNodeIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	return m_ptrName[uNodeIndex];
	}

unsigned Tree::GetLeafId(unsigned uNodeIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(IsLeaf(uNodeIndex));
	return m_Ids[uNodeIndex];
	}

// Append a new branch.
// This adds two new nodes and joins them to an existing leaf node.
// Return value is k, new nodes have indexes k and k+1 respectively.
unsigned Tree::AppendBranch(unsigned uExistingLeafIndex)
	{
	if (0 == m_uNodeCount)
		Quit("Tree::AppendBranch: tree has not been created");

#if	DEBUG
	assert(uExistingLeafIndex < m_uNodeCount);
	if (!IsLeaf(uExistingLeafIndex))
		{
		LogMe();
		Quit("AppendBranch(%u): not leaf", uExistingLeafIndex);
		}
#endif

	if (m_uNodeCount >= m_uCacheCount - 2)
		ExpandCache();

	const unsigned uNewLeaf1 = m_uNodeCount;
	const unsigned uNewLeaf2 = m_uNodeCount + 1;

	m_uNodeCount += 2;

	assert(m_uNeighbor2[uExistingLeafIndex] == NULL_NEIGHBOR);
	assert(m_uNeighbor3[uExistingLeafIndex] == NULL_NEIGHBOR);

	m_uNeighbor2[uExistingLeafIndex] = uNewLeaf1;
	m_uNeighbor3[uExistingLeafIndex] = uNewLeaf2;

	m_uNeighbor1[uNewLeaf1] = uExistingLeafIndex;
	m_uNeighbor1[uNewLeaf2] = uExistingLeafIndex;

	m_uNeighbor2[uNewLeaf1] = NULL_NEIGHBOR;
	m_uNeighbor2[uNewLeaf2] = NULL_NEIGHBOR;

	m_uNeighbor3[uNewLeaf1] = NULL_NEIGHBOR;
	m_uNeighbor3[uNewLeaf2] = NULL_NEIGHBOR;

	m_dEdgeLength2[uExistingLeafIndex] = 0;
	m_dEdgeLength3[uExistingLeafIndex] = 0;

	m_dEdgeLength1[uNewLeaf1] = 0;
	m_dEdgeLength2[uNewLeaf1] = 0;
	m_dEdgeLength3[uNewLeaf1] = 0;

	m_dEdgeLength1[uNewLeaf2] = 0;
	m_dEdgeLength2[uNewLeaf2] = 0;
	m_dEdgeLength3[uNewLeaf2] = 0;

	m_bHasEdgeLength1[uNewLeaf1] = false;
	m_bHasEdgeLength2[uNewLeaf1] = false;
	m_bHasEdgeLength3[uNewLeaf1] = false;

	m_bHasEdgeLength1[uNewLeaf2] = false;
	m_bHasEdgeLength2[uNewLeaf2] = false;
	m_bHasEdgeLength3[uNewLeaf2] = false;

	m_bHasHeight[uNewLeaf1] = false;
	m_bHasHeight[uNewLeaf2] = false;

	m_Ids[uNewLeaf1] = uInsane;
	m_Ids[uNewLeaf2] = uInsane;
	return uNewLeaf1;
	}

void Tree::LogMe() const
	{
	Log("Tree::LogMe %u nodes, ", m_uNodeCount);

	if (IsRooted())
		{
		Log("rooted.\n");
		Log("\n");
		Log("Index  Parnt  LengthP  Left   LengthL  Right  LengthR     Id  Name\n");
		Log("-----  -----  -------  ----   -------  -----  -------  -----  ----\n");
		}
	else
		{
		Log("unrooted.\n");
		Log("\n");
		Log("Index  Nbr_1  Length1  Nbr_2  Length2  Nbr_3  Length3     Id  Name\n");
		Log("-----  -----  -------  -----  -------  -----  -------  -----  ----\n");
		}

	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		Log("%5u  ", uNodeIndex);
		const unsigned n1 = m_uNeighbor1[uNodeIndex];
		const unsigned n2 = m_uNeighbor2[uNodeIndex];
		const unsigned n3 = m_uNeighbor3[uNodeIndex];
		if (NULL_NEIGHBOR != n1)
			{
			Log("%5u  ", n1);
			if (m_bHasEdgeLength1[uNodeIndex])
				Log("%7.4f  ", m_dEdgeLength1[uNodeIndex]);
			else
				Log("      *  ");
			}
		else
			Log("                ");

		if (NULL_NEIGHBOR != n2)
			{
			Log("%5u  ", n2);
			if (m_bHasEdgeLength2[uNodeIndex])
				Log("%7.4f  ", m_dEdgeLength2[uNodeIndex]);
			else
				Log("      *  ");
			}
		else
			Log("                ");

		if (NULL_NEIGHBOR != n3)
			{
			Log("%5u  ", n3);
			if (m_bHasEdgeLength3[uNodeIndex])
				Log("%7.4f  ", m_dEdgeLength3[uNodeIndex]);
			else
				Log("      *  ");
			}
		else
			Log("                ");

		if (m_Ids != 0 && IsLeaf(uNodeIndex))
			{
			unsigned uId = m_Ids[uNodeIndex];
			if (uId == uInsane)
				Log("    *");
			else
				Log("%5u", uId);
			}
		else
			Log("     ");

		if (m_bRooted && uNodeIndex == m_uRootNodeIndex)
			Log("  [ROOT] ");
		const char *ptrName = m_ptrName[uNodeIndex];
		if (ptrName != 0)
			Log("  %s", ptrName);
		Log("\n");
		}
	}

void Tree::SetEdgeLength(unsigned uNodeIndex1, unsigned uNodeIndex2,
  double dLength)
	{
	assert(uNodeIndex1 < m_uNodeCount && uNodeIndex2 < m_uNodeCount);
	assert(IsEdge(uNodeIndex1, uNodeIndex2));

	if (m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
		{
		m_dEdgeLength1[uNodeIndex1] = dLength;
		m_bHasEdgeLength1[uNodeIndex1] = true;
		}
	else if (m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
		{
		m_dEdgeLength2[uNodeIndex1] = dLength;
		m_bHasEdgeLength2[uNodeIndex1] = true;

		}
	else
		{
		assert(m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
		m_dEdgeLength3[uNodeIndex1] = dLength;
		m_bHasEdgeLength3[uNodeIndex1] = true;
		}

	if (m_uNeighbor1[uNodeIndex2] == uNodeIndex1)
		{
		m_dEdgeLength1[uNodeIndex2] = dLength;
		m_bHasEdgeLength1[uNodeIndex2] = true;
		}
	else if (m_uNeighbor2[uNodeIndex2] == uNodeIndex1)
		{
		m_dEdgeLength2[uNodeIndex2] = dLength;
		m_bHasEdgeLength2[uNodeIndex2] = true;
		}
	else
		{
		assert(m_uNeighbor3[uNodeIndex2] == uNodeIndex1);
		m_dEdgeLength3[uNodeIndex2] = dLength;
		m_bHasEdgeLength3[uNodeIndex2] = true;
		}
	}

unsigned Tree::UnrootFromFile()
	{
#if	TRACE
	Log("Before unroot:\n");
	LogMe();
#endif

	if (!m_bRooted)
		Quit("Tree::Unroot, not rooted");

// Convention: root node is always node zero
	assert(IsRoot(0));
	assert(NULL_NEIGHBOR == m_uNeighbor1[0]);

	const unsigned uThirdNode = m_uNodeCount++;

	m_uNeighbor1[0] = uThirdNode;
	m_uNeighbor1[uThirdNode] = 0;

	m_uNeighbor2[uThirdNode] = NULL_NEIGHBOR;
	m_uNeighbor3[uThirdNode] = NULL_NEIGHBOR;

	m_dEdgeLength1[0] = 0;
	m_dEdgeLength1[uThirdNode] = 0;
	m_bHasEdgeLength1[uThirdNode] = true;

	m_bRooted = false;

#if	TRACE
	Log("After unroot:\n");
	LogMe();
#endif

	return uThirdNode;
	}

// In an unrooted tree, equivalent of GetLeft/Right is 
// GetFirst/SecondNeighbor.
// uNeighborIndex must be a known neighbor of uNodeIndex.
// This is the way to find the other two neighbor nodes of
// an internal node.
// The labeling as "First" and "Second" neighbor is arbitrary.
// Calling these functions on a leaf returns NULL_NEIGHBOR, as
// for GetLeft/Right.
unsigned Tree::GetFirstNeighbor(unsigned uNodeIndex, unsigned uNeighborIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(uNeighborIndex < m_uNodeCount);
	assert(IsEdge(uNodeIndex, uNeighborIndex));

	for (unsigned n = 0; n < 3; ++n)
		{
		unsigned uNeighbor = GetNeighbor(uNodeIndex, n);
		if (NULL_NEIGHBOR != uNeighbor && uNeighborIndex != uNeighbor)
			return uNeighbor;
		}
	return NULL_NEIGHBOR;
	}

unsigned Tree::GetSecondNeighbor(unsigned uNodeIndex, unsigned uNeighborIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(uNeighborIndex < m_uNodeCount);
	assert(IsEdge(uNodeIndex, uNeighborIndex));

	bool bFoundOne = false;
	for (unsigned n = 0; n < 3; ++n)
		{
		unsigned uNeighbor = GetNeighbor(uNodeIndex, n);
		if (NULL_NEIGHBOR != uNeighbor && uNeighborIndex != uNeighbor)
			{
			if (bFoundOne)
				return uNeighbor;
			else
				bFoundOne = true;
			}
		}
	return NULL_NEIGHBOR;
	}

// Compute the number of leaves in the sub-tree defined by an edge
// in an unrooted tree. Conceptually, the tree is cut at this edge,
// and uNodeIndex2 considered the root of the sub-tree.
unsigned Tree::GetLeafCountUnrooted(unsigned uNodeIndex1, unsigned uNodeIndex2,
  double *ptrdTotalDistance) const
	{
	assert(!IsRooted());

	if (IsLeaf(uNodeIndex2))
		{
		*ptrdTotalDistance = GetEdgeLength(uNodeIndex1, uNodeIndex2);
		return 1;
		}

// Recurse down the rooted sub-tree defined by cutting the edge
// and considering uNodeIndex2 as the root.
	const unsigned uLeft = GetFirstNeighbor(uNodeIndex2, uNodeIndex1);
	const unsigned uRight = GetSecondNeighbor(uNodeIndex2, uNodeIndex1);

	double dLeftDistance;
	double dRightDistance;

	const unsigned uLeftCount = GetLeafCountUnrooted(uNodeIndex2, uLeft,
	  &dLeftDistance);
	const unsigned uRightCount = GetLeafCountUnrooted(uNodeIndex2, uRight,
	  &dRightDistance);

	*ptrdTotalDistance = dLeftDistance + dRightDistance;
	return uLeftCount + uRightCount;
	}

void Tree::RootUnrootedTree(ROOT Method)
	{
	assert(!IsRooted());
#if TRACE
	Log("Tree::RootUnrootedTree, before:");
	LogMe();
#endif

	unsigned uNode1;
	unsigned uNode2;
	double dLength1;
	double dLength2;
	FindRoot(*this, &uNode1, &uNode2, &dLength1, &dLength2, Method);

	if (m_uNodeCount == m_uCacheCount)
		ExpandCache();
	m_uRootNodeIndex = m_uNodeCount++;

	double dEdgeLength = GetEdgeLength(uNode1, uNode2);

	m_uNeighbor1[m_uRootNodeIndex] = NULL_NEIGHBOR;
	m_uNeighbor2[m_uRootNodeIndex] = uNode1;
	m_uNeighbor3[m_uRootNodeIndex] = uNode2;

	if (m_uNeighbor1[uNode1] == uNode2)
		m_uNeighbor1[uNode1] = m_uRootNodeIndex;
	else if (m_uNeighbor2[uNode1] == uNode2)
		m_uNeighbor2[uNode1] = m_uRootNodeIndex;
	else
		{
		assert(m_uNeighbor3[uNode1] == uNode2);
		m_uNeighbor3[uNode1] = m_uRootNodeIndex;
		}

	if (m_uNeighbor1[uNode2] == uNode1)
		m_uNeighbor1[uNode2] = m_uRootNodeIndex;
	else if (m_uNeighbor2[uNode2] == uNode1)
		m_uNeighbor2[uNode2] = m_uRootNodeIndex;
	else
		{
		assert(m_uNeighbor3[uNode2] == uNode1);
		m_uNeighbor3[uNode2] = m_uRootNodeIndex;
		}

	OrientParent(uNode1, m_uRootNodeIndex);
	OrientParent(uNode2, m_uRootNodeIndex);

	SetEdgeLength(m_uRootNodeIndex, uNode1, dLength1);
	SetEdgeLength(m_uRootNodeIndex, uNode2, dLength2);

	m_bHasHeight[m_uRootNodeIndex] = false;

	m_ptrName[m_uRootNodeIndex] = 0;

	m_bRooted = true;

#if	TRACE
	Log("\nPhy::RootUnrootedTree, after:");
	LogMe();
#endif

	Validate();
	}

bool Tree::HasEdgeLength(unsigned uNodeIndex1, unsigned uNodeIndex2) const
	{
	assert(uNodeIndex1 < m_uNodeCount);
	assert(uNodeIndex2 < m_uNodeCount);
	assert(IsEdge(uNodeIndex1, uNodeIndex2));

	if (m_uNeighbor1[uNodeIndex1] == uNodeIndex2)
		return m_bHasEdgeLength1[uNodeIndex1];
	else if (m_uNeighbor2[uNodeIndex1] == uNodeIndex2)
		return m_bHasEdgeLength2[uNodeIndex1];
	assert(m_uNeighbor3[uNodeIndex1] == uNodeIndex2);
	return m_bHasEdgeLength3[uNodeIndex1];
	}

void Tree::OrientParent(unsigned uNodeIndex, unsigned uParentNodeIndex)
	{
	if (NULL_NEIGHBOR == uNodeIndex)
		return;

	if (m_uNeighbor1[uNodeIndex] == uParentNodeIndex)
		;
	else if (m_uNeighbor2[uNodeIndex] == uParentNodeIndex)
		{
		double dEdgeLength2 = m_dEdgeLength2[uNodeIndex];
		m_uNeighbor2[uNodeIndex] = m_uNeighbor1[uNodeIndex];
		m_dEdgeLength2[uNodeIndex] = m_dEdgeLength1[uNodeIndex];
		m_uNeighbor1[uNodeIndex] = uParentNodeIndex;
		m_dEdgeLength1[uNodeIndex] = dEdgeLength2;
		}
	else
		{
		assert(m_uNeighbor3[uNodeIndex] == uParentNodeIndex);
		double dEdgeLength3 = m_dEdgeLength3[uNodeIndex];
		m_uNeighbor3[uNodeIndex] = m_uNeighbor1[uNodeIndex];
		m_dEdgeLength3[uNodeIndex] = m_dEdgeLength1[uNodeIndex];
		m_uNeighbor1[uNodeIndex] = uParentNodeIndex;
		m_dEdgeLength1[uNodeIndex] = dEdgeLength3;
		}

	OrientParent(m_uNeighbor2[uNodeIndex], uNodeIndex);
	OrientParent(m_uNeighbor3[uNodeIndex], uNodeIndex);
	}

unsigned Tree::FirstDepthFirstNode() const
	{
	assert(IsRooted());

// Descend via left branches until we hit a leaf
	unsigned uNodeIndex = m_uRootNodeIndex;
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetLeft(uNodeIndex);
	return uNodeIndex;
	}

unsigned Tree::FirstDepthFirstNodeR() const
	{
	assert(IsRooted());

// Descend via left branches until we hit a leaf
	unsigned uNodeIndex = m_uRootNodeIndex;
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetRight(uNodeIndex);
	return uNodeIndex;
	}

unsigned Tree::NextDepthFirstNode(unsigned uNodeIndex) const
	{
#if	TRACE
	Log("NextDepthFirstNode(%3u) ", uNodeIndex);
#endif

	assert(IsRooted());
	assert(uNodeIndex < m_uNodeCount);

	if (IsRoot(uNodeIndex))
		{
#if	TRACE
		Log(">> Node %u is root, end of traversal\n", uNodeIndex);
#endif
		return NULL_NEIGHBOR;
		}

	unsigned uParent = GetParent(uNodeIndex);
	if (GetRight(uParent) == uNodeIndex)
		{
#if	TRACE
		Log(">> Is right branch, return parent=%u\n", uParent);
#endif
		return uParent;
		}

	uNodeIndex = GetRight(uParent);
#if	TRACE
		Log(">> Descend left from right sibling=%u ... ", uNodeIndex);
#endif
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetLeft(uNodeIndex);

#if	TRACE
	Log("bottom out at leaf=%u\n", uNodeIndex);
#endif
	return uNodeIndex;
	}

unsigned Tree::NextDepthFirstNodeR(unsigned uNodeIndex) const
	{
#if	TRACE
	Log("NextDepthFirstNode(%3u) ", uNodeIndex);
#endif

	assert(IsRooted());
	assert(uNodeIndex < m_uNodeCount);

	if (IsRoot(uNodeIndex))
		{
#if	TRACE
		Log(">> Node %u is root, end of traversal\n", uNodeIndex);
#endif
		return NULL_NEIGHBOR;
		}

	unsigned uParent = GetParent(uNodeIndex);
	if (GetLeft(uParent) == uNodeIndex)
		{
#if	TRACE
		Log(">> Is left branch, return parent=%u\n", uParent);
#endif
		return uParent;
		}

	uNodeIndex = GetLeft(uParent);
#if	TRACE
		Log(">> Descend right from left sibling=%u ... ", uNodeIndex);
#endif
	while (!IsLeaf(uNodeIndex))
		uNodeIndex = GetRight(uNodeIndex);

#if	TRACE
	Log("bottom out at leaf=%u\n", uNodeIndex);
#endif
	return uNodeIndex;
	}

void Tree::UnrootByDeletingRoot()
	{
	assert(IsRooted());
	assert(m_uNodeCount >= 3);

	const unsigned uLeft = GetLeft(m_uRootNodeIndex);
	const unsigned uRight = GetRight(m_uRootNodeIndex);

	m_uNeighbor1[uLeft] = uRight;
	m_uNeighbor1[uRight] = uLeft;

	bool bHasEdgeLength = HasEdgeLength(m_uRootNodeIndex, uLeft) &&
	  HasEdgeLength(m_uRootNodeIndex, uRight);
	if (bHasEdgeLength)
		{
		double dEdgeLength = GetEdgeLength(m_uRootNodeIndex, uLeft) +
		  GetEdgeLength(m_uRootNodeIndex, uRight);
		m_dEdgeLength1[uLeft] = dEdgeLength;
		m_dEdgeLength1[uRight] = dEdgeLength;
		}

// Remove root node entry from arrays
	const unsigned uMoveCount = m_uNodeCount - m_uRootNodeIndex;
	const unsigned uUnsBytes = uMoveCount*sizeof(unsigned);
	memmove(m_uNeighbor1 + m_uRootNodeIndex, m_uNeighbor1 + m_uRootNodeIndex + 1,
	  uUnsBytes);
	memmove(m_uNeighbor2 + m_uRootNodeIndex, m_uNeighbor2 + m_uRootNodeIndex + 1,
	  uUnsBytes);
	memmove(m_uNeighbor3 + m_uRootNodeIndex, m_uNeighbor3 + m_uRootNodeIndex + 1,
	  uUnsBytes);

	const unsigned uDoubleBytes = uMoveCount*sizeof(double);
	memmove(m_dEdgeLength1 + m_uRootNodeIndex, m_dEdgeLength1 + m_uRootNodeIndex + 1,
	  uDoubleBytes);
	memmove(m_dEdgeLength2 + m_uRootNodeIndex, m_dEdgeLength2 + m_uRootNodeIndex + 1,
	  uDoubleBytes);
	memmove(m_dEdgeLength3 + m_uRootNodeIndex, m_dEdgeLength3 + m_uRootNodeIndex + 1,
	  uDoubleBytes);

	const unsigned uBoolBytes = uMoveCount*sizeof(bool);
	memmove(m_bHasEdgeLength1 + m_uRootNodeIndex, m_bHasEdgeLength1 + m_uRootNodeIndex + 1,
	  uBoolBytes);
	memmove(m_bHasEdgeLength2 + m_uRootNodeIndex, m_bHasEdgeLength2 + m_uRootNodeIndex + 1,
	  uBoolBytes);
	memmove(m_bHasEdgeLength3 + m_uRootNodeIndex, m_bHasEdgeLength3 + m_uRootNodeIndex + 1,
	  uBoolBytes);

	const unsigned uPtrBytes = uMoveCount*sizeof(char *);
	memmove(m_ptrName + m_uRootNodeIndex, m_ptrName + m_uRootNodeIndex + 1, uPtrBytes);

	--m_uNodeCount;
	m_bRooted = false;

// Fix up table entries
	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
#define DEC(x)	if (x != NULL_NEIGHBOR && x > m_uRootNodeIndex) --x;
		DEC(m_uNeighbor1[uNodeIndex])
		DEC(m_uNeighbor2[uNodeIndex])
		DEC(m_uNeighbor3[uNodeIndex])
#undef	DEC
		}

	Validate();
	}

unsigned Tree::GetLeafParent(unsigned uNodeIndex) const
	{
	assert(IsLeaf(uNodeIndex));

	if (IsRooted())
		return GetParent(uNodeIndex);

	if (m_uNeighbor1[uNodeIndex] != NULL_NEIGHBOR)
		return m_uNeighbor1[uNodeIndex];
	if (m_uNeighbor2[uNodeIndex] != NULL_NEIGHBOR)
		return m_uNeighbor2[uNodeIndex];
	return m_uNeighbor3[uNodeIndex];
	}

// TODO: This is not efficient for large trees, should cache.
double Tree::GetNodeHeight(unsigned uNodeIndex) const
	{
	if (!IsRooted())
		Quit("Tree::GetNodeHeight: undefined unless rooted tree");
	
	if (IsLeaf(uNodeIndex))
		return 0.0;

	if (m_bHasHeight[uNodeIndex])
		return m_dHeight[uNodeIndex];

	const unsigned uLeft = GetLeft(uNodeIndex);
	const unsigned uRight = GetRight(uNodeIndex);
	double dLeftLength = GetEdgeLength(uNodeIndex, uLeft);
	double dRightLength = GetEdgeLength(uNodeIndex, uRight);

	if (dLeftLength < 0)
		dLeftLength = 0;
	if (dRightLength < 0)
		dRightLength = 0;

	const double dLeftHeight = dLeftLength + GetNodeHeight(uLeft);
	const double dRightHeight = dRightLength + GetNodeHeight(uRight);
	const double dHeight = (dLeftHeight + dRightHeight)/2;
	m_bHasHeight[uNodeIndex] = true;
	m_dHeight[uNodeIndex] = dHeight;
	return dHeight;
	}

unsigned Tree::GetNeighborSubscript(unsigned uNodeIndex, unsigned uNeighborIndex) const
	{
	assert(uNodeIndex < m_uNodeCount);
	assert(uNeighborIndex < m_uNodeCount);
	if (uNeighborIndex == m_uNeighbor1[uNodeIndex])
		return 0;
	if (uNeighborIndex == m_uNeighbor2[uNodeIndex])
		return 1;
	if (uNeighborIndex == m_uNeighbor3[uNodeIndex])
		return 2;
	return NULL_NEIGHBOR;
	}

unsigned Tree::GetNeighbor(unsigned uNodeIndex, unsigned uNeighborSubscript) const
	{
	switch (uNeighborSubscript)
		{
	case 0:
		return m_uNeighbor1[uNodeIndex];
	case 1:
		return m_uNeighbor2[uNodeIndex];
	case 2:
		return m_uNeighbor3[uNodeIndex];
		}
	Quit("Tree::GetNeighbor, sub=%u", uNeighborSubscript);
	return NULL_NEIGHBOR;
	}

// TODO: check if this is a performance issue, could cache a lookup table
unsigned Tree::LeafIndexToNodeIndex(unsigned uLeafIndex) const
	{
	const unsigned uNodeCount = GetNodeCount();
	unsigned uLeafCount = 0;
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (IsLeaf(uNodeIndex))
			{
			if (uLeafCount == uLeafIndex)
				return uNodeIndex;
			else
				++uLeafCount;
			}
		}
	Quit("LeafIndexToNodeIndex: out of range");
	return 0;
	}

unsigned Tree::GetLeafNodeIndex(const char *ptrName) const
	{
	const unsigned uNodeCount = GetNodeCount();
	for (unsigned uNodeIndex = 0; uNodeIndex < uNodeCount; ++uNodeIndex)
		{
		if (!IsLeaf(uNodeIndex))
			continue;
		const char *ptrLeafName = GetLeafName(uNodeIndex);
		if (0 == strcmp(ptrName, ptrLeafName))
			return uNodeIndex;
		}
	Quit("Tree::GetLeafNodeIndex, name not found");
	return 0;
	}

void Tree::Copy(const Tree &tree)
	{
	const unsigned uNodeCount = tree.GetNodeCount();
	InitCache(uNodeCount);

	m_uNodeCount = uNodeCount;

	const size_t UnsignedBytes = uNodeCount*sizeof(unsigned);
	const size_t DoubleBytes = uNodeCount*sizeof(double);
	const size_t BoolBytes = uNodeCount*sizeof(bool);

	memcpy(m_uNeighbor1, tree.m_uNeighbor1, UnsignedBytes);
	memcpy(m_uNeighbor2, tree.m_uNeighbor2, UnsignedBytes);
	memcpy(m_uNeighbor3, tree.m_uNeighbor3, UnsignedBytes);

	memcpy(m_Ids, tree.m_Ids, UnsignedBytes);

	memcpy(m_dEdgeLength1, tree.m_dEdgeLength1, DoubleBytes);
	memcpy(m_dEdgeLength2, tree.m_dEdgeLength2, DoubleBytes);
	memcpy(m_dEdgeLength3, tree.m_dEdgeLength3, DoubleBytes);

	memcpy(m_dHeight, tree.m_dHeight, DoubleBytes);

	memcpy(m_bHasEdgeLength1, tree.m_bHasEdgeLength1, BoolBytes);
	memcpy(m_bHasEdgeLength2, tree.m_bHasEdgeLength2, BoolBytes);
	memcpy(m_bHasEdgeLength3, tree.m_bHasEdgeLength3, BoolBytes);

	memcpy(m_bHasHeight, tree.m_bHasHeight, BoolBytes);

	m_uRootNodeIndex = tree.m_uRootNodeIndex;
	m_bRooted = tree.m_bRooted;

	for (unsigned uNodeIndex = 0; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		if (tree.IsLeaf(uNodeIndex))
			{
			const char *ptrName = tree.GetLeafName(uNodeIndex);
			m_ptrName[uNodeIndex] = strsave(ptrName);
			}
		else
			m_ptrName[uNodeIndex] = 0;
		}

#if	DEBUG
	Validate();
#endif
	}

// Create rooted tree from a vector description.
// Node indexes are 0..N-1 for leaves, N..2N-2 for
// internal nodes.
// Vector subscripts are i-N and have values for
// internal nodes only, but those values are node
// indexes 0..2N-2. So e.g. if N=6 and Left[2]=1,
// this means that the third internal node (node index 8)
// has the second leaf (node index 1) as its left child.
// uRoot gives the vector subscript of the root, so add N
// to get the node index.
void Tree::Create(unsigned uLeafCount, unsigned uRoot, const unsigned Left[],
  const unsigned Right[], const float LeftLength[], const float RightLength[],
  const unsigned LeafIds[], char **LeafNames)
	{
	Clear();

	m_uNodeCount = 2*uLeafCount - 1;
	InitCache(m_uNodeCount);

	for (unsigned uNodeIndex = 0; uNodeIndex < uLeafCount; ++uNodeIndex)
		{
		m_Ids[uNodeIndex] = LeafIds[uNodeIndex];
		m_ptrName[uNodeIndex] = strsave(LeafNames[uNodeIndex]);
		}

	for (unsigned uNodeIndex = uLeafCount; uNodeIndex < m_uNodeCount; ++uNodeIndex)
		{
		unsigned v = uNodeIndex - uLeafCount;
		unsigned uLeft = Left[v];
		unsigned uRight = Right[v];
		float fLeft = LeftLength[v];
		float fRight = RightLength[v];

		m_uNeighbor2[uNodeIndex] = uLeft;
		m_uNeighbor3[uNodeIndex] = uRight;

		m_bHasEdgeLength2[uNodeIndex] = true;
		m_bHasEdgeLength3[uNodeIndex] = true;

		m_dEdgeLength2[uNodeIndex] = fLeft;
		m_dEdgeLength3[uNodeIndex] = fRight;

		m_uNeighbor1[uLeft] = uNodeIndex;
		m_uNeighbor1[uRight] = uNodeIndex;

		m_dEdgeLength1[uLeft] = fLeft;
		m_dEdgeLength1[uRight] = fRight;

		m_bHasEdgeLength1[uLeft] = true;
		m_bHasEdgeLength1[uRight] = true;
		}

	m_bRooted = true;
	m_uRootNodeIndex = uRoot + uLeafCount;

	Validate();
	}
