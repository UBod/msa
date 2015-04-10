#include "muscle.h"
#include "clust.h"

void Clust::InsertMetric(unsigned uIndex1, unsigned uIndex2, float dMetric)
	{
	RBInsert(uIndex1, uIndex2, dMetric);
	}

void Clust::DeleteMetric(unsigned uIndex)
	{
	for (unsigned uNodeIndex = GetFirstCluster(); uNodeIndex != uInsane;
	  uNodeIndex = GetNextCluster(uNodeIndex))
		{
		if (uIndex == uNodeIndex)
			continue;
		DeleteMetric(uIndex, uNodeIndex);
		}
	}

void Clust::InitMetric(unsigned uMaxNodeIndex)
	{
	m_uRBNodeCount = m_uTriangularMatrixSize;
	m_RBParent = new unsigned[m_uRBNodeCount];
	m_RBLeft = new unsigned[m_uRBNodeCount];
	m_RBRight = new unsigned[m_uRBNodeCount];
	m_RBi = new ushort[m_uRBNodeCount];
	m_RBj = new ushort[m_uRBNodeCount];
	m_RBMetric = new float[m_uRBNodeCount];
	m_RBColor = new bool[m_uRBNodeCount];
	m_RBRoot = RB_NIL;

#if	DEBUG
	{
// Initialize fields to invalid values so we have a chance
// catch attempts to use them if they're not properly set.
	unsigned InvalidNode = m_uRBNodeCount + 1;
	for (unsigned Node = 0; Node < m_uRBNodeCount; ++Node)
		{
		m_RBParent[Node] = InvalidNode;
		m_RBLeft[Node] = InvalidNode;
		m_RBRight[Node] = InvalidNode;
		m_RBi[Node] = InvalidNode;
		m_RBj[Node] = InvalidNode;
		}
	}
#endif
	}

void Clust::ListMetric() const
	{
	Log("Red-black tree root=%u\n", m_RBRoot);
	Log("\n");
	Log(" Node  Parent   Left  Right  Color      i      j  Metric\n");
	Log("-----  ------  -----  -----  -----  -----  -----  ------\n");

	if (RB_NIL == m_RBRoot)
		return;

	unsigned Count = 0;
	unsigned Start = RBMin(m_RBRoot);
	for (unsigned Node = Start; RB_NIL != Node; Node = RBNext(Node))
		{
		Log("%5u", Node);

		if (RB_NIL != m_RBParent[Node])
			Log("  %6u", m_RBParent[Node]);
		else
			Log("        ");

		if (RB_NIL != m_RBLeft[Node])
			Log("  %5u", m_RBLeft[Node]);
		else
			Log("       ");

		if (RB_NIL != m_RBRight[Node])
			Log("  %5u", m_RBRight[Node]);
		else
			Log("       ");

		Log("  %s  %5u  %5u  %g\n",
		  m_RBColor[Node] ? "  Red" : "Black",
		  m_RBi[Node],
		  m_RBj[Node],
		  m_RBMetric[Node]);

		if (++Count > m_uRBNodeCount)
			{
			Log(" ** LOOP ** \n");
			break;
			}
		}
	}

// If there is a left subtree, predecessor is the
// largest key found under the left branch. Otherwise,
// is first node in path to root that is a right child.
unsigned Clust::RBPrev(unsigned Node) const
	{
	assert(Node < m_uRBNodeCount);

	unsigned Left = m_RBLeft[Node];
	if (RB_NIL != Left)
		return RBMax(Left);

	for (;;)
		{
		unsigned Parent = m_RBParent[Node];
		if (RB_NIL == Parent)
			return RB_NIL;
		if (m_RBRight[Parent] == Node)
			return Parent;
		Node = Parent;
		}
	}

// If there is a right subtree, sucessor is the
// smallest key found under the right branch. Otherwise,
// is first node in path to root that is a left child.
unsigned Clust::RBNext(unsigned Node) const
	{
	if (Node >= m_uRBNodeCount)
		Quit("RBNext(%u)", Node);
	assert(Node < m_uRBNodeCount);

	unsigned Right = m_RBRight[Node];
	if (RB_NIL != Right)
		return RBMin(Right);

	for (;;)
		{
		unsigned Parent = m_RBParent[Node];
		if (RB_NIL == Parent)
			return RB_NIL;
		if (m_RBLeft[Parent] == Node)
			return Parent;
		Node = Parent;
		}
	}

// Minimum is in leftmost leaf
unsigned Clust::RBMin(unsigned RBNode) const
	{
	assert(RB_NIL != RBNode);
	for (;;)
		{
		unsigned Left = m_RBLeft[RBNode];
		if (RB_NIL == Left)
			return RBNode;
		RBNode = Left;
		}
	}

// Maximum is in rightmost leaf
unsigned Clust::RBMax(unsigned RBNode) const
	{
	assert(RB_NIL != RBNode);
	for (;;)
		{
		unsigned Right = m_RBRight[RBNode];
		if (RB_NIL == Right)
			return RBNode;
		RBNode = Right;
		}
	}

void Clust::DeleteMetric(unsigned uIndex1, unsigned uIndex2)
	{
	unsigned RBNode = (unsigned) VectorIndex(uIndex1, uIndex2);
	RBDelete(RBNode);
	}

void Clust::RBDelete(unsigned Node)
	{
#if	DEBUG
	ValidateRB();
	//Log("@@ Before RBDelete(%u)\n", Node);
	//ListMetric();
#endif

	unsigned Left = m_RBLeft[Node];
	unsigned Right = m_RBRight[Node];
	unsigned Parent = m_RBParent[Node];

// If one or two nil children, splice out this node.
	if (RB_NIL == Left || RB_NIL == Right)
		{
//		Log("@@ One child\n");
	// Child is non-NIL child, or NIL if none.
		unsigned Child = (Left != RB_NIL ? Left : Right);

	// Special case if root
		if (RB_NIL == Parent)
			{
			assert(Node == m_RBRoot);
			m_RBRoot = Child;
			if (RB_NIL != Child)
				m_RBParent[Child] = RB_NIL;
			return;
			}

	// Typical case.
	// Update parent->child link
		if (m_RBLeft[Parent] == Node)
			m_RBLeft[Parent] = Child;
		else
			{
			assert(m_RBRight[Parent] == Node);
			m_RBRight[Parent] = Child;
			}

	// Update child->parent link
		if (RB_NIL != Child)
			m_RBParent[Child] = Parent;

#if	DEBUG
		//Log("@@ After RBDelete(%u)\n", Node);
		//ListMetric();
		ValidateRB();
#endif
		return;
		}

	//Log("@@ RBDelete(%u) Tricky case\n", Node);
	//ListMetric();

// Trickier case, node has two children.
	assert(Left != RB_NIL && Right != RB_NIL);

// We're going to splice out successor node from its
// current position and insert it in place of node
// to be deleted.

// Successor cannot be nil because there is a right child.
	unsigned Next = RBNext(Node);
	assert(Next != RB_NIL);

// The successor of a node with two children is
// guaranteed to have no more than one child.
	unsigned NextLeft = m_RBLeft[Next];
	unsigned NextRight = m_RBRight[Next];
	assert(RB_NIL == NextLeft || RB_NIL == NextRight);

// Successor of node with two children cannot be the root.
	unsigned NextParent = m_RBParent[Next];
	assert(RB_NIL != NextParent);

// Ugly special case if successor is right child
	if (Next == Right)
		{
#if DEBUG
		//Log("@@ Before RBDelete(%u) (tricky next==right)\n", Node);
		//ListMetric();
#endif
		m_RBParent[Next] = Parent;

		if (RB_NIL == Parent)
			{
			m_RBRoot = Next;
			m_RBParent[Next] = RB_NIL;
			}
		else
			{
			if (m_RBLeft[Parent] == Node)
				m_RBLeft[Parent] = Next;
			else
				{
				assert(m_RBRight[Parent] == Node);
				m_RBRight[Parent] = Next;
				}
			}

		m_RBLeft[Next] = Left;

		if (RB_NIL != Left)
			m_RBParent[Left] = Next;

#if	DEBUG
		//Log("@@ After RBDelete(%u) (tricky next==right)\n", Node);
		//ListMetric();
		ValidateRB();
#endif
		return;
		}

// Set NextChild either to the one child of successor, or nil.
	unsigned NextChild = (NextLeft != RB_NIL ? NextLeft : NextRight);

// Splice successor from its current position
	if (m_RBLeft[NextParent] == Next)
		m_RBLeft[NextParent] = NextChild;
	else
		{
		assert(m_RBRight[NextParent] == Next);
		m_RBRight[NextParent] = NextChild;
		}

	if (RB_NIL != NextChild)
		m_RBParent[NextChild] = NextParent;

// Insert successor into position currently held by node
// to be deleted.
	if (RB_NIL == Parent)
		{
		m_RBRoot = Next;
		m_RBParent[Next] = RB_NIL;
		}
	else
		{
		if (m_RBLeft[Parent] == Node)
			m_RBLeft[Parent] = Next;
		else
			{
			assert(m_RBRight[Parent] == Node);
			m_RBRight[Parent] = Next;
			}
		}

	m_RBLeft[Next] = Left;
	m_RBRight[Next] = Right;
	m_RBParent[Next] = Parent;

	m_RBParent[Left] = Next;
	m_RBParent[Right] = Next;

#if	DEBUG
	//Log("@@ After RBDelete(%u)\n", Node);
	//ListMetric();
	ValidateRB();
#endif
	}

unsigned Clust::RBInsert(unsigned i, unsigned j, float fMetric)
	{
#if	DEBUG
	ValidateRB();
#endif

	unsigned NewNode = VectorIndex(i, j);
	m_RBMetric[NewNode] = fMetric;
	m_RBi[NewNode] = i;
	m_RBj[NewNode] = j;

// New node is always inserted as a leaf.
// Proof that this is possible is found in algorithm
// textbooks (I forget the argument).
	m_RBLeft[NewNode] = RB_NIL;
	m_RBRight[NewNode] = RB_NIL;

	unsigned NewParent = RB_NIL;
	unsigned Node = m_RBRoot;

	unsigned uCount = 0;
	while (RB_NIL != Node)
		{
		NewParent = Node;
		if (fMetric < m_RBMetric[Node])
			Node = m_RBLeft[Node];
		else
			Node = m_RBRight[Node];
		++uCount;
		if (uCount > m_uRBNodeCount)
			Quit("Infinite loop in RBInsert");
		}

	m_RBParent[NewNode] = NewParent;
	if (RB_NIL == NewParent)
		m_RBRoot = NewNode;
	else
		{
		if (fMetric < m_RBMetric[NewParent])
			m_RBLeft[NewParent] = NewNode;
		else
			m_RBRight[NewParent] = NewNode;
		}

#if	DEBUG
	{
	unsigned Next = RBNext(NewNode);
	if (Next != RB_NIL)
		assert(NewNode == RBPrev(Next));
	unsigned Prev = RBPrev(NewNode);
	if (Prev != RB_NIL)
		assert(NewNode == RBNext(Prev));
	ValidateRB();
	}
#endif
	return NewNode;
	}

void Clust::ValidateRBNode(unsigned Node, const char szMsg[]) const
	{
	if (RB_NIL == Node)
		return;

	unsigned Parent = m_RBParent[Node];
	unsigned Left = m_RBLeft[Node];
	unsigned Right = m_RBRight[Node];

	unsigned Next = RBNext(Node);
	unsigned Prev = RBPrev(Node);

	if (RB_NIL != Next && RBPrev(Next) != Node)
		{
		ListMetric();
		Quit("ValidateRB(%s) Node=%u Next=%u Prev(Next)=%u",
		  szMsg, Node, Next, RBPrev(Next));
		}

	if (RB_NIL != Prev && RBNext(Prev) != Node)
		{
		ListMetric();
		Quit("ValidateRB(%s) Node=%u Prev=%u Next(Prev)=%u",
		  szMsg, Node, Prev, RBNext(Prev));
		}

	if (RB_NIL != Parent)
		{
		if (m_RBLeft[Parent] != Node && m_RBRight[Parent] != Node)
			{
			ListMetric();
			Quit("ValidateRB(%s): Parent %u not linked to child %u\n",
			  szMsg, Parent, Node);
			}
		}

	if (RB_NIL != Left)
		{
		if (m_RBParent[Left] != Node)
			{
			ListMetric();
			Quit("ValidateRB(%s): Left child %u not linked to parent %u\n",
			  szMsg, Left, Node);
			}
		}

	if (RB_NIL != Right)
		{
		if (m_RBParent[Right] != Node)
			{
			ListMetric();
			Quit("ValidateRB(%s): Right child %u not linked to parent %u\n",
			  szMsg, Right, Node);
			}
		}

	ValidateRBNode(Left, szMsg);
	ValidateRBNode(Right, szMsg);
	}

void Clust::ValidateRB(const char szMsg[]) const
	{
	if (RB_NIL == m_RBRoot)
		return;

	ValidateRBNode(m_RBRoot, szMsg);

	unsigned Node = RBMin(m_RBRoot);
	for (;;)
		{
		unsigned Next = RBNext(Node);
		if (RB_NIL == Next)
			break;
		if (m_RBMetric[Node] > m_RBMetric[Next])
			{
			ListMetric();
			Quit("ValidateRBNode(%s): metric out of order %u=%g %u=%g",
			  szMsg, Node, m_RBMetric[Node], Next, m_RBMetric[Next]);
			}
		Node = Next;
		}
	}
