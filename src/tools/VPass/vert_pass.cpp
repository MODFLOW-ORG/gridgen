
#include "vert_pass.h"
//int vpnum = 0;

void determineArea(VPassBox& vpb)
{
	double area = (double)(FLT_MAX);
	double r = fabs(vpb.bSrc->dx * vpb.bSrc->dy) ;
	if(r < area)
		area  =r;

	r = fabs(vpb.bDst->dx * vpb.bDst->dy);
	if(r< area)
		area = r;
	for(list<Box*>::iterator bit = vpb.bPssList.begin(); bit != vpb.bPssList.end(); ++bit)
	{
		Box* tb = *bit;
		r = fabs(tb->dx * tb->dy);
		if( r< area)
			area = r;
	}
	vpb.area = area;
}
void genOneVP(list<Box*>& bl, vector<VPassBox>& vpbl)
{
	if(bl.size() >= 3)
	{
		//assert(bl.front()->active && bl.back()->active);
		VPassBox vpb;
		vpb.bSrc = bl.front();
		bl.pop_front();

		vpb.bDst = bl.back();
		bl.pop_back();

		vpb.bPssList.insert(vpb.bPssList.end(), bl.begin(), bl.end());
		determineArea(vpb);

		vpbl.push_back(vpb);

		//vpnum++;
		//cout<<vpnum<<"----vertical pass--id :"<<vpb.bSrc->id<<"\t"<<vpb.bDst->id<<"\t"<<vpb.bPssList.front()->id<<endl;
	}
}
//exploration in two vertical directions
void verticalExplore(Box* tb, list<Box*>& bl, vector<VPassBox>& vpbl,  bool topToDown = true)
{
	list<Box*> nei;

	if(topToDown && tb->pDown != NULL)
		tb->pDown->getLeaves(nei);
	else if( (!topToDown) && (tb->pUp != NULL))
		tb->pUp->getLeaves(nei);


	for(list<Box*>::iterator bit = nei.begin(); bit != nei.end();++bit)
	{
		Box* pb = *bit;
		assert(pb != tb);
		if(pb == NULL || !(pb->isLeaf))
			continue;
		if(pb->active)
		{
			if(bl.size() >= 2)
			{
				list<Box*> tlist(bl.begin(), bl.end());
				tlist.push_back(pb);
				genOneVP(tlist, vpbl);
			}
		}
		else
		{
			list<Box*> tlist(bl.begin(), bl.end());
			tlist.push_back(pb);
			verticalExplore(pb, tlist, vpbl, topToDown);
		}
	}
}

void boxVPConnection(Box* b, vector<VPassBox>& vpbl, bool topToDown)
{
	if(b == NULL || (!(b->isLeaf)) || (b->active == false))
		return;

	list<Box*> bl;
	bl.push_back(b);
	verticalExplore(b, bl, vpbl,  topToDown);
}

void getVPT(Box* b, vector<VPTPair>& vpts, bool topToDown)
{
	if(b == NULL || (!(b->isLeaf)) || (b->active == false))
		return;

	vector<VPassBox> vpbs;
	boxVPConnection(b, vpbs, topToDown);
	for(vector<VPassBox>::iterator vit = vpbs.begin(); vit != vpbs.end(); ++vit)
	{
		Box* b1 = vit->bSrc; 
		Box* b2 = vit->bDst;

		vector<VPassBox>::iterator jit = vit;
		++jit;
		for(; jit != vpbs.end(); )
		{
			if(b1 == jit->bSrc && b2 == jit->bDst)
			{
				vit->area += jit->area;
				jit = vpbs.erase(jit);
			}
			else
			{
				++jit;
			}
		}
	}

	//coy the final results
	for(vector<VPassBox>::iterator vit = vpbs.begin(); vit != vpbs.end(); ++vit)
	{
		VPTPair vpt(vit->bSrc, vit->bDst, vit->area);
		vpts.push_back(vpt);
	}

	////////////////////////////////////////////
	//for testing
	//cout<<"total vertical pass through : "<<vpnum<<endl;
}

////this function 
//void VerticalPass(Grid* grid, vector<VPassBox>& vpbl)
//{
//	//cast
//	QuadTree3D* qtgrid=dynamic_cast<QuadTree3D*>(grid);
//	ModflowGrid* mfgrid=dynamic_cast<ModflowGrid*>(grid);
//	ModflowGrid2D* mf2dqtgrid=dynamic_cast<ModflowGrid2D*>(grid);
//
//	if(qtgrid == NULL && mfgrid == NULL)
//		return;
//
//	//vector<VPassBox> vpbl;
//	Box* tb = NULL;
//	for(NodeGroup::iterator nit = grid->nodegroup.begin(); nit != grid->nodegroup.end(); ++nit)
//	{
//		tb = *nit;	
//
//		if(tb == NULL || (!(tb->isLeaf)) || (tb->active == false))
//			continue;
//
//		boxVPConnection(tb, vpbl, false, true);
//	}
//
//	cerr<<"total vertical pass number is "<<vpbl.size()<<endl;
//
//	for(vector<VPassBox>::iterator vit = vpbl.begin(); vit !=vpbl.end(); ++vit)
//	{
//		vit->bSrc->flag = 1;
//
//		if(vit->bSrc->id == 22911)
//		{
//			cout<<"debugging here!\n";
//		}
//
//		vector<VPassBox> vbs;
//		boxVPConnection(vit->bSrc, vbs, false, true);
//		cout<<"id : "<<vit->bSrc->id<<"\t"<<vit->bPssList.front()->id<<"\t"<<vit->bDst->id<<endl;
//		//cout<<vit->bSrc->id<<"\t";
//	}
//	cout<<endl;
//	
//}
