#include "BoxPlyClip.h"
#include "intersection.h"
#include <map>
#include <fstream>
using namespace std;


namespace cusg
{
	bool compareTwoV(VNode* v1, VNode* v2)
	{
		if(fabs(v1->x - v2->x) >SMALLCOORDINATE)
		{
			cerr<<"coordiante not equal"<<endl;
			return false;
		}
		if(fabs(v1->y- v2->y) >SMALLCOORDINATE)
		{
			cerr<<"coordiante not equal"<<endl;
			return false;
		}
		if(v1->isInDirect != v2->isInDirect || v1->isOnline != v2->isOnline || v1->isScanPnt != v2->isScanPnt || v1->isVisited != v2->isVisited)
		{
			cerr<<"coordinate not equal"<<endl;
			return false;
		}
			return true;
	}
	bool compareTwoPly(VNode* h1, VNode* h2)
	{
		VNode* v1=h1,*v2=h2;
		do
		{
			if(!compareTwoV(v1, v2))
				return false;
			v1=v1->next; v2=v2->next;
		} while (v1 != h1 );
		if(v2 != h2)
			return false;
		return true;
	}
	bool compareTwoBPC(BoxPlyClip* b1, BoxPlyClip* b2)
	{
		list<VNode*>::iterator vit1 = b1->plys.begin();
		list<VNode*>::iterator vit2 = b2->plys.begin();
		for (; vit1 != b1->plys.end(); ++vit1, ++vit2)
		{
			if(!compareTwoPly(*vit1, *vit2))
				return false;
		}
		if(vit2 != b2->plys.end())
			cout<<"the ply size is not equal"<<endl;
		return true;
	}
	Point2d getCenter(VNode* head)
	{
		double x = 0, y=0;
		int size = 0;
		VNode* tmv = head;
		do{
			x += tmv->x;
			y+= tmv->y;
			size++;

			tmv = tmv->next;
		}while(tmv != head);

		return Point2d(x/size, y/size);
	}

	bool AlmostSame(double a[2], double b[2])
	{
		return (fabs(a[0] - b[0]) <= SMALLCOORDINATE && fabs(a[1] - b[1]) <= SMALLCOORDINATE);
	}
	bool AlmostEqual(VNode* v1, VNode* v2)
	{
		return (fabs(v1->x - v2->x) <=SMALLCOORDINATE && fabs(v1->y - v2->y) <= SMALLCOORDINATE);
	}
	//void removeDuplicate(VNode* vhead)
	//{
	//	VNode* vn = vhead;
	//	do 
	//	{
	//		if(vn->next!= vhead&&AlmostEqual(vn, vn->next))
	//		{
	//			VNode* tmp = vn->next;
	//			VNode* vnn = tmp->next;
	//			vn->next = vnn;
	//			vnn->prev = vn;
	//			delete tmp;
	//		}
	//		vn = vn->next;
	//	} while (vn != vhead);
	//}
	//void getFirstUnEqual(VNode* sv, 
	//sort the vertices based on X or Y ascendingly
	void InsertionSortV(vector<VNode*>& vs,bool isX)
	{
		int vsize = vs.size();
		if(vsize < 2)
			return;

		for (int i = 1; i < vsize; i++)
		{
			VNode* tv = vs[i];
			int j = i - 1;
			bool bcon = (!isX && (vs[j]->x > tv->x) ) ||(isX && (vs[j]->y > tv->y) ) ;
			while(j >=0 && ( (!isX && (vs[j]->x > tv->x) ) ||(isX && (vs[j]->y > tv->y) )))
			{
				vs[j+1] = vs[j];
				j--;
			}
			vs[j + 1] = tv;
		}
	}

	int locateVInVector(vector<VNode*>& vs, VNode* vv, int startIdx, int totalSize)
	{
		for (int i = startIdx; i < totalSize; i++)
		{
			if(vs[i] == vv)
				return i;
		}
		return -1;
	}

	void getXYBound(VNode* first, double& xmin, double& xmax, double& ymin, double& ymax)
	{
		xmin = ymin = FLT_MAX;
		xmax = ymax = FLT_MIN;
		VNode* vn = first;
		do 
		{
			if(vn->x < xmin)		xmin = vn->x;
			if(vn->y <ymin)		ymin = vn->y;
			if(vn->x > xmax)	xmax = vn->x;
			if(vn->y > ymax)	ymax = vn->y;
			vn = vn->next;
		} while (vn != first);
	}

	VNode* copyOnePly(VNode* head, bool bSetVisited)
	{
		if(head == NULL)
			return NULL;

		VNode* nh = new VNode(head->x, head->y);
		if(bSetVisited)
			head->isVisited = true;
		nh->next = nh->prev = nh;
		VNode* cv = head->next;
		while(cv != head)
		{
			VNode* t = new VNode(cv->x, cv->y);
			nh->prev->next = t;
			t->prev = nh->prev;
			t->next = nh;
			nh->prev = t;

			if(bSetVisited)
				cv->isVisited = true;

			cv = cv->next;
		}

		return nh;
	}
	void resetBPCStatus(VNode* head)
	{
		if(head == NULL)
			return ;
		VNode* vn = head;
		do 
		{
			vn->resetStatus();
			vn = vn->next;
		} while (vn != head);
	}
	BoxPlyClip* copyBoxPlyClip(BoxPlyClip* bpc, bool bSetVisited = false)
	{
		BoxPlyClip* nbcp = new BoxPlyClip();
		
		VNode* nh = copyOnePly(bpc->first, bSetVisited);
		nbcp->first = nh;
		nbcp->plys.push_back(nh);

		if(bpc->plys.size() > 1)
		{
			list<VNode*>::iterator vit =bpc->plys.begin();
			++vit;
			for(;vit != bpc->plys.end(); ++vit)
			{
				VNode* nvv = copyOnePly(*vit, bSetVisited);
				nbcp->plys.push_back(nvv);
			}
		}
		return nbcp;
	}
	void BoxPlyClip::resetStatus()
	{
		//delete the new added
		for(list<VNode*>::iterator vit = newAdds.begin(); vit != newAdds.end(); ++vit)
		{
			VNode* vn = *vit;
			VNode* pre1 = vn->prev;
			VNode* next1= vn->next;
			pre1->next = next1;
			next1->prev =pre1;
			vn->prev = vn->next = NULL;
			delete vn;
			//vit = newAdds.erase(vit);
		}
		newAdds.clear();

		//for(list<VNode*>::iterator vit = plys.begin(); vit != plys.end(); ++vit)
		//{
		//	resetBPCStatus(*vit);
		//}
		//loop over the influenced vertices
		for(list<VNode*>::iterator vit = influVs.begin(); vit != influVs.end(); ++vit)
		{
			VNode* vn = *vit;
			vn->resetStatus();
		}
		influVs.clear();
	}

	double getRotAngle(Vector2d& v1, Vector2d& v2)
	{
		//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
		//have bug here
		double angle = 0;
		Vector2d vc1 = v1.normalize();
		Vector2d vc2 = v2.normalize();
		double dotpro = (vc1[0] * vc2[0]) + (vc1[1] * vc2[1]);
		if(fabs(dotpro - 1.0) <= EPSION)
			return 0;
		else if(fabs(dotpro + 1.0) <= EPSION)
			return PI;
		else
		{
			double cross = 0;
			angle = acos(dotpro);
			cross = (vc1[0] * vc2[1]) - (vc2[0] * vc1[1]);
			if(cross < 0)
				angle = 2* PI - angle;
			
			return angle;
		}
	}
	//test enclosing
	bool ClipPlyEnclose(VNode* head, const Point2d& p2d)
	{
		double wnum = 0.0;
		if(head==NULL || head->next == head)
			return false;
		VNode* cv = head;
		do 
		{
			Vector2d vc1(cv->x - p2d[0], cv->y - p2d[1]);
			Vector2d vc2(cv->next->x - p2d[0], cv->next->y - p2d[1]);

			double ang = getRotAngle(vc1, vc2);
			if(ang > PI)
				ang = ang - 2* PI;
			wnum += ang;

			cv = cv->next;
		} while (cv != head);
		
		//if( (abs(wnum) < EPSION) || (abs(wnum - 2*PI) < EPSION) || (abs(wnum + 2*PI) < EPSION))
		if(fabs(wnum - 2* PI) < EPSION || fabs(wnum + 2* PI) <EPSION)	
			return true;

		return false;
	}

	void getFeasiblePreNext(VNode* tv, VNode* & vn1, VNode*& vn2, VNode*& vnprev1, VNode*& vnprev2)
	{
		vn1 = tv->prev;
		 vn2 = tv->next;
		while(AlmostEqual(vn1, tv))
		{
			vn1 = vn1->prev;
		}
		while(AlmostEqual(vn2, tv))
		{
			vn2 = vn2->next;
		}
		vnprev1 = vn1;
		vnprev2 = vn2;
		//filtering 
		while(vn1->isScanPnt)
		{
			vnprev1 = vn1;
			vn1 = vn1->prev;
		}
		while(vn2->isScanPnt)
		{
			vnprev2 = vn2;
			vn2 = vn2->next;
		}

	}
	bool decideNextOrderForOut(VNode* tv, bool isX, bool isBigger)
	{
		assert(tv->isInDirect);
		assert(tv->next->isInDirect || tv->prev->isInDirect);

		if(!(tv->prev->isInDirect && tv->next->isInDirect))
			return tv->next->isInDirect;
		else
		{
			//1:pre; 2:next
			VNode* vn1=NULL,*vn2=NULL,*vnpre1=NULL,*vnpre2=NULL;
			getFeasiblePreNext(tv, vn1, vn2, vnpre1, vnpre2);
			assert(!(vnpre1->isVisited&&vnpre2->isVisited));
			//if one of them is not in indirection
			if(!(vn1->isInDirect && vn2->isInDirect))
			{
				if( vn1->isInDirect && (!vnpre1->isVisited))//previous
					return false;
				else if(vn2->isInDirect && (!vnpre2->isVisited))//next
					return true;
				else
					assert(false);
			}
			//if both are in direction, but one of them has been visited
			if(vn1->isInDirect && vn2->isInDirect)
			{
				if(vnpre1->isVisited || vn1->isVisited)
					return true;
				if(vnpre2->isVisited || vn2->isVisited)
					return false;
			}
			//if(vnpre1->isVisited)
			//	return true;
			//else
			//	return false;
			if(!(vnpre1->isScanPnt && vnpre2->isScanPnt))
			{
				Vector2d vc1(vnpre1->x - tv->x, vnpre1->y - tv->y);
				Vector2d vc2(vnpre2->x - tv->x, vnpre2->y - tv->y);
				double angle = getRotAngle(vc1, vc2);
				if( (!isX && isBigger) || (isX && !isBigger) )
					return (angle <= PI);
				else 
					return (angle >= PI);
			}
			else
			{
				Vector2d vc1(vnpre1->x - tv->x, vnpre1->y - tv->y);
				Vector2d vc2(vn2->x - tv->x, vn2->y - tv->y);
				double angle = getRotAngle(vc1, vc2);
				if( (!isX && isBigger) || (isX && !isBigger) )
					return (angle <= PI);
				else 
					return (angle >= PI);
			}
		}

		//comment 10/4/2013
		//Vector2d vc1(tv->prev->x - tv->x, tv->prev->y - tv->y);
		//Vector2d vc2(tv->next->x - tv->x, tv->next->y - tv->y);
		//double angle = getRotAngle(vc1, vc2);
		//if( (!isX && isBigger) || (isX && !isBigger) )
		//	return (angle <= PI);
		//else 
		//	return (angle >= PI);
	}
	bool decideNextOrderForHole(VNode* tv, bool isX, bool isBigger)
	{
		assert(tv->isInDirect);
		assert(tv->next->isInDirect || tv->prev->isInDirect);

		if(!(tv->prev->isInDirect && tv->next->isInDirect))
			return tv->next->isInDirect;
		else
		{
			//1:pre; 2:next
			VNode* vn1=NULL,*vn2=NULL,*vnpre1=NULL,*vnpre2=NULL;
			getFeasiblePreNext(tv, vn1, vn2, vnpre1, vnpre2);
			assert(!(vnpre1->isVisited&&vnpre2->isVisited));
			if(!(vn1->isInDirect && vn2->isInDirect))
			{
				if( vn1->isInDirect && (!vnpre1->isVisited))
					return false;
				else if(vn2->isInDirect && (!vnpre2->isVisited))
					return true;
				else
					assert(false);
			}
			//if both are in direction, but one of them has been visited
			if(vn1->isInDirect && vn2->isInDirect)
			{
				if(vnpre1->isVisited || vn1->isVisited)
					return true;
				if(vnpre2->isVisited || vn2->isVisited)
					return false;
			}

			//if(vnpre1->isVisited)
			//	return true;
			//else
			//	return false;
		}

		Vector2d vc1(tv->prev->x - tv->x, tv->prev->y - tv->y);
		Vector2d vc2(tv->next->x - tv->x, tv->next->y - tv->y);

		double angle = getRotAngle(vc1, vc2);
		if( ( (!isX) && isBigger ) || (isX && (!isBigger) ) )
			return (angle >= PI);
		else
			return (angle <= PI);
		

	}

	bool checkAllIndirect(double fval, bool isX, bool isBigger,VNode* head)
	{

		Point2d p2d = getCenter(head);
		if(isX)
		{
			return (isBigger ==(p2d[0] > fval));
		}
		else
		{
			return (isBigger ==(p2d[1] > fval));
		}
	}
	BoxPlyClip::BoxPlyClip()
	{
		first = NULL;
	}
	BoxPlyClip::~BoxPlyClip()
	{
		destroy();
	}
	void BoxPlyClip::destroy()
	{
		for(list<VNode*>::iterator vit = plys.begin(); vit != plys.end(); ++vit)
		{
			destroy(*vit);
		}
		plys.clear();
		first = NULL;
	}
	void BoxPlyClip::destroy(VNode* vn)
	{
		if(vn==NULL)
			return;
		VNode* tmp = vn->next;
		while(tmp != vn)
		{
			VNode* tf = tmp;
			tmp = tmp->next;
			if(tf != NULL)
				delete tf;
		}
		delete vn;
		vn = NULL;
	}

	void BoxPlyClip::add(VNode* v)
	{
		if(first != NULL&& fabs(v->x - first->x) < SMALLNUMBER && fabs(v->y - first->y) < SMALLNUMBER)
			return;
		if(first != NULL&&fabs(v->x - first->prev->x) < SMALLNUMBER && fabs(v->y - first->prev->y) < SMALLNUMBER)
			return;

			if(first == NULL){
				first = v;
				first->next = v;
				first->prev = v;
			}
			else{
				VNode* next = first;
				VNode* prev = next->prev;
				next->prev = v;
				v->next = next;
				v->prev = prev;
				prev->next = v;
			}
	}
	void BoxPlyClip::insert(VNode* s, VNode* e, VNode* v)
	{
		//if(s != NULL&& abs(v->x - s->x) < SMALLNUMBER && abs(v->y - s->y) < SMALLNUMBER)
		//	return;
		//if(e != NULL&&abs(v->x - e->prev->x) < SMALLNUMBER && abs(v->y - e->prev->y) < SMALLNUMBER)
		//	return;

		s->next = v;
		v->prev = s;
		v->next = e;
		e->prev = v;
	}
	VNode* BoxPlyClip::collectAlongWay(vector<VNode*>& vs, VNode* tmpv, bool nextO)
	{
		assert(tmpv->isInDirect && tmpv->isScanPnt);
		VNode* ptv = tmpv;
		vs.push_back(ptv);
		VNode* cv = nextO ? tmpv->next : tmpv->prev;
		//ptv =cv;

		while(AlmostEqual(ptv, cv) && cv->isInDirect)
		{
			vs.push_back(cv);
			ptv = cv;
			cv = nextO ? cv->next : cv->prev;
		}
		//collect all the online vertices
		//while(cv->isScanPnt)
		//{
		//	vs.push_back(cv);
		//	cv = nextO ? cv->next : cv->prev;
		//}
	//	if(cv->isInDirect)
		///{
			while(cv->isInDirect)
			{
				vs.push_back(cv);
				if(cv->isScanPnt)
					break;
				cv = nextO ? cv->next : cv->prev;
			}
		//}

		return vs.back();

	}

	//set intersection
	int BoxPlyClip::setIntesection(VNode* startV, double fval, bool isX, bool isBigger, vector<VNode*>& intersectV)
	{
		double xmin = 0,ymin = 0, xmax = 0, ymax = 0;
		getXYBound(startV, xmin, xmax, ymin, ymax);
		double xwidth = fabs(xmax - xmin) ;
		double ywidth = fabs(ymax - ymin);
		if(xwidth < SMALLNUMBER)
			xwidth = ywidth;
		if(ywidth < SMALLNUMBER) 
			ywidth = xwidth;

		double a[2]={0},b[2]={0};
		if(isX)
		{
			a[0] = fval; a[1] = ymin - ywidth / 2;
			b[0] = fval; b[1]= ymax +ywidth/ 2;
		}
		else
		{
			a[0] = xmin- xwidth/2; a[1] = fval;
			b[0] = xmax+xwidth/2; b[1]= fval;
		}
		
		//try to intersect
		map<VNode*, bool> shouldRecord;
		vector<VNode*> interCand;
		VNode* stmp = startV;
		do 
		{
			 double c[2]={0}, d[2]={0}, p[2]={FLT_MIN, FLT_MIN};
			stmp->setDirection(fval, isX, isBigger);
			//intersect
			c[0] = stmp->x; c[1] = stmp->y;
			d[0] = stmp->next->x; d[1]= stmp->next->y;
			char code = SegSegInt(a, b, c, d, p);

			//set influenced vertices
			if(stmp->isInDirect)
				influVs.push_back(stmp);
			//if(!stmp->isOnline)
			//{
				if(code =='e')
				{
					stmp->isOnline = stmp->next->isOnline =true;
					stmp->isScanPnt = stmp->next->isScanPnt = true;//
					influVs.push_back(stmp);	influVs.push_back(stmp->next);
				}
				else if(code =='v')
				{
					//if(AlmostEqual(c, p))
					if(AlmostSame(c, p))
					{
						stmp->isOnline = true;	stmp->isScanPnt = true;//
						influVs.push_back(stmp);
					}
					else{
						//assert(AlmostEqual(d, p));
						assert(AlmostSame(d, p));
						stmp->next->isOnline = true; stmp->next->isScanPnt = true;//
						influVs.push_back(stmp->next);
					}
				}
				else if(code =='1')
				{
					//if(AlmostEqual(c, p))
					if(AlmostSame(c, p))
					{
						stmp->isOnline = true;	stmp->isScanPnt = true;//
						influVs.push_back(stmp);
					}
					//else if(AlmostEqual(d, p))
					else if(AlmostSame(d, p))
					{
						stmp->next->isOnline = true; stmp->next->isScanPnt = true;//
						influVs.push_back(stmp->next);
					}
					else
					{
						VNode* nv = new VNode(p[0], p[1]);
						this->insert(stmp, stmp->next, nv);
						//add to the newAdds array
						newAdds.push_back(nv);

						nv->isInDirect = true;
						nv->isScanPnt = true;//

						//interCand.push_back(nv);
						shouldRecord.insert(make_pair(nv, true));
						stmp = stmp->next;
					}
				}
			//}

			stmp = stmp->next;
		} while (stmp != startV);

		//loop to decide online v
		stmp = startV;
		do 
		{
			if(stmp->isOnline)
			{
				stmp->isInDirect = true;
				//stmp->isScanPnt = true;
				interCand.push_back(stmp);
				shouldRecord.insert(make_pair(stmp, true));
			}
			stmp = stmp->next;
		} while (stmp != startV);

		vector<VNode*> uselessNodes;
		vector<VNode*> notOnline;
		//vector<VNode*> uselessNodes;
		for(vector<VNode*>::iterator vit =interCand.begin(); vit != interCand.end(); ++vit)
		{
			VNode* vn = *vit;
			
			bool pureInDir1 = (vn->next->isInDirect && (!vn->next->isOnline));
			bool pureInDir2 = (vn->prev->isInDirect &&(!vn->prev->isOnline));
			if(vn->isOnline && (!(pureInDir1&&pureInDir2)) )
			{
				notOnline.push_back(vn);
				//vn->isOnline = false;
			}

			if(pureInDir1 || pureInDir2)
				shouldRecord[vn] = true;
			else
				shouldRecord[vn] = false;

			if(!(pureInDir1 || pureInDir2))
			{
				uselessNodes.push_back(vn);
				//vn->isInDirect =false; vn->isScanPnt = false;
			}
		}
		//not online
		for(vector<VNode*>::iterator vit = notOnline.begin(); vit != notOnline.end(); ++vit)
		{
			VNode* vn = *vit;
			vn->isOnline = false;
		}
		//useless nodes
		for(vector<VNode*>::iterator vit = uselessNodes.begin(); vit != uselessNodes.end(); ++vit)
		{
			VNode* vn = *vit;
			vn->isOnline = false;
			vn->isInDirect = false;// vn->isScanPnt = false;
		}


		//decide the final intersection vertices
		int onNum = 0;
		for(map<VNode*, bool>::iterator mit = shouldRecord.begin(); mit != shouldRecord.end(); ++mit)
		{
			bool b = mit->second;
			VNode* v = mit->first;
			if(b)
			{
				intersectV.push_back(v);
				onNum++;
			}
		}

		return onNum;
	}

	void BoxPlyClip::collectSecondPart(vector<VNode*>& interVs, int sIdx, int eIdx,VNode* sv, VNode* ev, bool isNextOrder, int vsize, list<BoxPlyClip*>& bpcResults,
		bool isX, bool isBigger, vector<VNode*>& secondHalf)
	{
		int startId = sIdx; 
		bool nextO = isNextOrder;

		//collect these vertices
		vector<VNode*> tmpVec;
		VNode* tsvTM = collectAlongWay(tmpVec, interVs[startId], nextO);

		for(vector<VNode*>::iterator vit = tmpVec.begin(); vit != tmpVec.end(); ++vit)
		{
			VNode* vn = *vit;
			if(!vn->isVisited)
			{
				secondHalf.push_back(vn); vn->isVisited = true;
			}
		}
		
		int tsvEndIdx = locateVInVector(interVs, tsvTM, startId, vsize);
		assert(tsvEndIdx != -1);

		//comment 10/4/2013
		//if(tsvEndIdx >=eIdx)
		//	return;

		/**************************************************************/
		//deal with the first part
		int fsIdx = startId + 1;
		bool fnextO = false;//decideNextOrderForOut(interVs[fsIdx], isX, isBigger);
		VNode* vsn1 = nextO ?  interVs[startId]->prev : interVs[startId]->next;
		if(interVs[startId]->isOnline && vsn1->isInDirect &&(!vsn1->isVisited))
		{
			fsIdx = startId; fnextO = !nextO;
		}
		if(fsIdx != startId && fsIdx < tsvEndIdx)
			fnextO = decideNextOrderForOut(interVs[fsIdx], isX, isBigger);

		int feIdx =tsvEndIdx - 1;
		VNode* ven2 = nextO ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
		if(interVs[tsvEndIdx]->isOnline && ven2->isInDirect &&(!ven2->isVisited))
		{
				bool eo = decideNextOrderForHole(interVs[tsvEndIdx], isX, isBigger);// !(decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger));
				VNode* enext = eo ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
				if(enext->isVisited)
				{
					feIdx  = tsvEndIdx;
				}
		}

		//deal with the first part
		if(fsIdx <feIdx)
		{
			collectPolygonsInInterval(isX, isBigger, interVs, bpcResults, fsIdx, feIdx, fnextO);
		}
		//end of first part
		/**********************************************************************************************/


		/*********************************************************************************************/
		//deal with the second part
		int ssIdx = tsvEndIdx + 1;
		bool sorder = false;//decideNextOrderForOut(interVs[ssIdx], isX, isBigger);
		VNode* vn1 = nextO ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
		if(interVs[tsvEndIdx]->isOnline && vn1->isInDirect && (!vn1->isVisited))
		{
			bool eo = decideNextOrderForHole(interVs[tsvEndIdx], isX, isBigger); //!(decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger));
			//bool eo = decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger);
			VNode* snext = eo ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
			if(!snext->isVisited)
			{
				ssIdx =tsvEndIdx ; sorder = nextO;
			}
		}

		if(ssIdx != tsvEndIdx && ssIdx <eIdx)
		{
			sorder = decideNextOrderForOut(interVs[ssIdx], isX, isBigger);
		}

		//if(interVs[tsvEndIdx]->isOnline)
		//{
		//	bool eo = decideNextOrderForHole(interVs[tsvEndIdx], isX, isBigger); //!(decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger));
		//	//bool eo = decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger);
		//	VNode* snext = eo ? interVs[tsvEndIdx]->next : interVs[tsvEndIdx]->prev;
		//	if(!snext->isVisited)
		//	//if(snext->isVisited)
		//	{
		//		assert( interVs[tsvEndIdx]->next->isInDirect &&  interVs[tsvEndIdx]->prev->isInDirect);
		//		ssIdx =tsvEndIdx ; sorder = nextO;
		//	}
		//	//else
		//	//{
		//	//	bool eo = decideNextOrderForOut(interVs[tsvEndIdx], isX, isBigger);
		//	//	assert(false);
		//	//}
		//}
		if(ssIdx < eIdx)
		{
			collectSecondPart(interVs, ssIdx, eIdx, interVs[ssIdx], ev, sorder, vsize, bpcResults, isX, isBigger, secondHalf);
		}
		//end of second part
		/*********************************************************************************************/

	}
	
	void BoxPlyClip::constructClipPly(vector<VNode*>& firstHalf, vector<VNode*>& secondHalf, list<BoxPlyClip*>& bpcResults)
	{
		if(firstHalf.size() + secondHalf.size() <3)
			return ;
		else
		{
			BoxPlyClip* bpc = new BoxPlyClip();
			for (vector<VNode*>::iterator vit = firstHalf.begin(); vit != firstHalf.end(); ++vit)
			{
				bpc->add(copy(*vit));
			}
			for (vector<VNode*>::reverse_iterator rvit = secondHalf.rbegin(); rvit != secondHalf.rend(); ++rvit)
			{
				bpc->add(copy(*rvit));
			}
			bpc->plys.push_back(bpc->first);

			c_ply cp;
			convertClipToPly(bpc->first, cp);
			if(fabs(cp.getArea()) < SMALLAREA)//SMALLNUMBER
			{
				cp.destroy();
				bpc->destroy();
				return;
			}
			else
			{
				cp.destroy();
				bpcResults.push_back(bpc);
			}

		}
	}
	void BoxPlyClip::collectOnePly(vector<VNode*>& interVs, int sIdx, int eIdx, bool isNextOrder, list<BoxPlyClip*>& bpcResults,
		bool isX, bool isBigger)
	{
		int vsize = interVs.size();
		vector<VNode*> firstHalf;
		vector<VNode*> secondHalf;

		//collect first half: out boundary
		VNode* sv =interVs[sIdx];
		VNode* ev = interVs[eIdx];
		VNode* tv = sv;
		firstHalf.push_back(sv);
		sv->isVisited = true;
		do{
			if(isNextOrder)		tv = tv->next;
			else	tv = tv->prev;

			tv->isVisited = true;
			firstHalf.push_back(tv);
		} while (tv != ev);

		//filter the starting vertex
		int startId =sIdx + 1;
		int endId = eIdx - 1;
		bool nextO = !isNextOrder;
		VNode* vno = !isNextOrder ? interVs[sIdx]->next : interVs[sIdx]->prev;
		if(interVs[sIdx]->isOnline && (!vno->isVisited))
		{
			assert(vno->isInDirect);
			nextO = !isNextOrder;
			startId = sIdx;
		}
		else
		{
			startId = sIdx + 1;
			if(startId <endId)
				nextO = decideNextOrderForOut(interVs[startId], isX, isBigger);

		}

		//if(interVs[eIdx]->isOnline

		if(startId < endId)
		{
			collectSecondPart(interVs, startId, endId, interVs[startId], interVs[endId], nextO,vsize, bpcResults, isX, isBigger, secondHalf);
		}

		//combine two lists to get a new polygon
		constructClipPly(firstHalf, secondHalf, bpcResults);
	}

	void BoxPlyClip::collectPolygonsInInterval(bool isX, bool isBigger, vector<VNode*>& interVs, list<BoxPlyClip*>& cres, int from, int to,bool ino)
	{
		int size = interVs.size();
		bool isNextO = ino;
		//VNode* sv = interVs[from];
		
		int startId = from;
		while(startId < to)
		{
			VNode* sv = interVs[startId];
			assert(sv->next->isInDirect || sv->prev->isInDirect);

			vector<VNode*> tvec;
			VNode* ev = collectAlongWay(tvec, sv, isNextO);
			int endId = locateVInVector(interVs, ev, startId, size);

			//assert(endId != -1 && startId != endId);
			assert(endId != -1);
		
			collectOnePly(interVs, startId, endId, isNextO, cres, isX, isBigger);//note here or not, may have bugs here

			if(endId >=to)
			{
				assert(endId == to);
				break;
			}

			VNode* vn1 = isNextO ? interVs[endId]->next : interVs[endId]->prev;
			VNode* vn2=  isNextO ? interVs[endId]->prev : interVs[endId]->prev;
			if(interVs[endId]->isOnline && (!vn1->isVisited))
			{
				startId = endId;
			}
			else
			{
				startId = endId + 1;
				isNextO = decideNextOrderForOut(interVs[startId], isX, isBigger);
			}
		}

	}

	void BoxPlyClip::collectPolygons(double fval, bool isX, bool isBigger, vector<int>& voNums, vector<VNode*>& interVs, list<BoxPlyClip*>& cres)
	{
		int size = interVs.size();
		if(size <=1)
		{
			if(checkAllIndirect(fval, isX, isBigger,first))
			{
				cres.push_back(copyBoxPlyClip(this, true));
			}
			//set visited
			//for(list<VNode*>::iterator )
			return;
		}

		//sort the intersection vertices
		InsertionSortV(interVs, isX);

		bool isNextO = decideNextOrderForOut(interVs[0], isX, isBigger);//need modify here

		collectPolygonsInInterval(isX, isBigger, interVs, cres, 0, size - 1, isNextO);

	}


	void BoxPlyClip::clip(double fval, bool isX, bool isBigger, list<BoxPlyClip*>& cres)
	{
		vector<int> onlineNums;
		vector<VNode*> interVs;
		//set intersections for holes
		for (list<VNode*>::iterator vit = plys.begin(); vit != plys.end(); ++vit)
		{
			VNode* sv = *vit;
			int vo = setIntesection(sv, fval, isX, isBigger,interVs);
			onlineNums.push_back(vo);
		}

		//collect the out boundary
		collectPolygons(fval, isX, isBigger, onlineNums, interVs, cres);

		/**********************************************/
		if(plys.size()==0 || onlineNums.size() == 0)
			return;

		list<VNode*>::iterator pit = plys.begin();
		vector<int>::iterator iit =onlineNums.begin();
		//++pit; ++iit; //discard the out boundary
		for (; pit != plys.end(); ++pit, ++iit)
		{
			
			if(*iit == 0 && ((*pit)->isVisited == false)  && (*pit)->isInDirect)
			{
				if((*pit)->isScanPnt)
				{
					VNode* vn1=NULL, *vn2=NULL,*vnprev1=NULL,*vnprev2=NULL;
					getFeasiblePreNext(*pit, vn1, vn2, vnprev1, vnprev2);
					if(!vnprev1->isInDirect || !vnprev2->isInDirect)
						continue;
				}

				bool haveAdded = false;
				VNode* nv = *pit;
				c_ply nply(c_ply::PIN);
				convertClipToPly(nv, nply);
				const Point2d& ctr = nply.getCenter();
				for (list<BoxPlyClip*>::iterator bit = cres.begin(); bit != cres.end(); ++bit)
				{
					if(ClipPlyEnclose((*bit)->first, ctr)){
						(*bit)->plys.push_back(copyOnePly(nv));
						haveAdded = true;
						break;
					}
				}
				if(!haveAdded)
				{
					BoxPlyClip* bpc = convertPlyToClip(nply);
					bpc->plys.push_back(bpc->first);
					cres.push_back(bpc);
				}
				nply.destroy();//clear memory
			}

		}
	}
	void BoxPlyClip::clip(double xmin, double xmax, double ymin, double ymax, list<BoxPlyClip*>& cres)
	{
		//pre-processing
		list<BoxPlyClip*> res1, res2, res3; 
		this->clip(xmin, true, true, res1);

		for(list<BoxPlyClip*>::iterator lit  = res1.begin(); lit != res1.end(); ++lit)
		{
			BoxPlyClip* bpc = *lit;
			bpc->clip(xmax, true, false, res2);

			bpc->destroy();
			delete bpc;
		}

		for(list<BoxPlyClip*>::iterator lit  = res2.begin(); lit != res2.end(); ++lit)
		{
			BoxPlyClip* bpc = *lit;
			bpc->clip(ymin, false, true, res3);
			
			bpc->destroy();
			delete bpc;
		}

		for(list<BoxPlyClip*>::iterator lit  = res3.begin(); lit != res3.end(); ++lit)
		{
			BoxPlyClip* bpc = *lit;
			bpc->clip(ymax, false, false, cres);

			bpc->destroy();
			delete bpc;
		}
		
	}

	VNode* BoxPlyClip::copy(const VNode* v)
	{
		VNode* nv = new VNode(v->x, v->y);
		return nv;
	}

	VNode* convertPlyToVN(const c_ply& tply)
	{
		ply_vertex* head = tply.getHead();
		VNode* nh = new VNode(head->getPos()[0], head->getPos()[1]);
		nh->next = nh->prev = nh;
		ply_vertex* cv = head->getNext();
		while(cv != head)
		{
			VNode* t = new VNode(cv->getPos()[0], cv->getPos()[1]);
			nh->prev->next = t;
			t->prev = nh->prev;
			t->next = nh;
			nh->prev = t;

			cv = cv->getNext();
		}

		return nh;
		
	}

	//convert the ply to BoxPlyClip
	BoxPlyClip* convertPlyToClip(const c_ply& tply)
	{
		BoxPlyClip* bpc = new BoxPlyClip();
		ply_vertex* thead = tply.getHead();
		ply_vertex* tv = thead;
		do 
		{
			VNode* tnode = new VNode(tv->getPos()[0], tv->getPos()[1]);
			bpc->add(tnode);
			
			tv = tv->getNext();
		} while (tv!=thead);
		return bpc;
	}

	//convert
	void convertClipToPly(VNode* fv, c_ply& resPly)
	{
		resPly.beginPoly();
		VNode* nf = fv;
		do 
		{
			resPly.addVertex(nf->x, nf->y,true);
			nf = nf->next;
		} while (nf != fv);

		resPly.endPoly(false);
		//resPly.endPoly(true);
	}

	//convert
	BoxPlyClip* convertCPolygonToClip(const GIS_polygon& polygon)
	{
		BoxPlyClip* bpc = new BoxPlyClip();
		for(GIS_polygon::const_iterator cit = polygon.begin(); cit != polygon.end(); ++cit)
		{
			const c_ply& ccp = *cit;
			if(ccp.getType()== c_ply::POUT)
			{
				VNode* nvn = convertPlyToVN(ccp);
				bpc->first = nvn;
				bpc->plys.push_front(nvn);
			}
			else
			{
				VNode* nvn = convertPlyToVN(ccp);
				bpc->plys.push_back(nvn);
			}
		}
		return bpc;
	}
	//convert
	void convertClipToPolygon(BoxPlyClip* bpc, GIS_polygon* cpolygon)
	{
		if(bpc == NULL || bpc->first == NULL)
			return;
		//loop over 
		c_ply outply(c_ply::POUT);
		convertClipToPly(bpc->first, outply);
		cpolygon->push_back(outply);

		if(bpc->plys.size() ==1)
			return;
		list<VNode*>::iterator vit = bpc->plys.begin();
		++vit;
		for(; vit != bpc->plys.end(); ++vit)
		{
			c_ply tp(c_ply::PIN);
			convertClipToPly(*vit, tp);
			cpolygon->push_back(tp);
		}
		
		
	}

BoxPlyClip* convertSimplePolygonToClip(const c_polygon& polygon)
{
		BoxPlyClip* bpc = new BoxPlyClip();
		for(c_polygon::const_iterator cit = polygon.begin(); cit != polygon.end(); ++cit)
		{
			const c_ply& ccp = *cit;
			if(ccp.getType()== c_ply::POUT)
			{
				VNode* nvn = convertPlyToVN(ccp);
				bpc->first = nvn;
				bpc->plys.push_front(nvn);
			}
			else
			{
				VNode* nvn = convertPlyToVN(ccp);
				bpc->plys.push_back(nvn);
			}
		}
		return bpc;
}

VNode* convertBoxToClip(Point2d ld, Point2d rd, Point2d ru, Point2d lu)
{
	VNode* tn[4];
	tn[0] = new VNode(ld[0], ld[1]);
	tn[1] = new VNode(rd[0], rd[1]);
	tn[2] = new VNode(ru[0], ru[1]);
	tn[3] = new VNode(lu[0], lu[1]);

	VNode* nh = tn[0];
	nh->next = nh->prev = tn[0];

	for(int i = 1; i <= 3; i++)
	{
		nh->prev->next = tn[i];
		tn[i]->prev = nh->prev;
		tn[i]->next = nh;
		nh->prev = tn[i];
	}

	return nh;
}

double getOutBdArea(BoxPlyClip* bpc)
{
	c_ply outply(c_ply::POUT);
    convertClipToPly(bpc->first, outply);

	double area = outply.getArea();
	outply.destroy();
	return area;
}

//make sure the min and max are in correct order
double intersectTwoBox(Point2d ld, Point2d rd, Point2d ru, Point2d lu, 
	double xmin2, double xmax2, double ymin2, double ymax2)
{
	double area(0.0);

	//build a BoxClip
	BoxPlyClip bpc;
	VNode* nh = convertBoxToClip(ld, rd, ru, lu);
	bpc.first = nh;
	bpc.plys.push_back(nh);
	
	list<BoxPlyClip*> cres;
	bpc.clip(xmin2, xmax2, ymin2, ymax2, cres);

	for (list<BoxPlyClip*>::iterator cit = cres.begin(); cit != cres.end(); ++cit)
	{
		//GIS_polygon nplygon;
		//convertClipToPolygon(*cit, &nplygon);
		//area += nplygon.getArea();

		area += getOutBdArea(*cit);

		(*cit)->destroy();
		delete *cit;
	}

	bpc.destroy();
	return area;
}

//void testClip()
//{
//	Point2d pt[4];
//	pt[0].set(0, 1);
//	pt[1].set(1, 0);
//	pt[2].set(2, 1);
//	pt[3].set(1, 2);
//
//	double a1 = intersectTwoBox(pt[0], pt[1], pt[2], pt[3], -1, 3, -1, 3);
//	double a2 = intersectTwoBox(pt[0], pt[1], pt[2], pt[3], 0.5, 1.5, 0.5, 1.5);
//	double a3 = intersectTwoBox(pt[0], pt[1], pt[2], pt[3], 0.75, 1.25, 0.75, 1.25);
//	double a4 = intersectTwoBox(pt[0], pt[1], pt[2], pt[3], 0.5, 1.5, 0, 2);
//
//	assert(false);
//}


}
