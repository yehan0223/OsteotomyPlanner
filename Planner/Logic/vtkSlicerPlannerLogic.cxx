/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// Planner Logic includes
#include "vtkSlicerPlannerLogic.h"

// Slicer CLI includes
#include <qSlicerCoreApplication.h>
#include <qSlicerModuleManager.h>
#include "qSlicerAbstractCoreModule.h"
#include <qSlicerCLIModule.h>
#include <vtkSlicerCLIModuleLogic.h>

//Slicer IO
#include <QFileInfo>
#include <QDir>
#include <QDateTime>
#include <ctkMessageBox.h>
#include <ctkUtils.h>
#include <QDebug>
#include "qMRMLUtils.h"
#include <vtkMRMLSceneViewNode.h>


// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLHierarchyNode.h>
#include <vtkMRMLModelHierarchyNode.h>
#include <vtkMRMLModelStorageNode.h>
#include <vtkMRMLModelDisplayNode.h>

// VTK includes
#include <vtkNew.h>
#include <vtkMassProperties.h>
#include <vtkTriangleFilter.h>
#include <vtkAppendPolyData.h>
#include "vtkVector.h"
#include "vtkVectorOperators.h"
#include "vtkMath.h"
#include "vtkCutter.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyDataNormals.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataPointSampler.h"
#include "vtkDecimatePro.h"
#include "vtkMatrix4x4.h"
#include "vtkVertexGlyphFilter.h"
#include <vtksys/SystemTools.hxx>

// STD includes
#include <cassert>
#include <sstream>


//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerPlannerLogic);

//----------------------------------------------------------------------------
//Constructor
vtkSlicerPlannerLogic::vtkSlicerPlannerLogic()
{
  this->SkullWrappedPreOP = NULL;
  this->HealthyBrain = NULL;
  this->BoneTemplate = NULL;
  this->splitLogic = NULL;
  this->wrapperLogic = NULL;
  this->preOPICV = 0;
  this->healthyBrainICV = 0;
  this->currentICV = 0;
  this->templateICV = 0;
  this->TempMerged = NULL;
  this->TempWrapped = NULL;
  this->CurrentModel = NULL;
  this->SourcePoints = NULL;
  this->SourcePointsDense = NULL;
  this->TargetPoints = NULL;
  this->Fiducials = NULL;
  this->cellLocator = NULL;
  this->bendMode = Double;
  this->bendSide = A;
  this->BendingPlane = NULL;
  this->BendingPlaneLocator = NULL;
  this->bendInitialized = false;
  this->BendingPolyData = NULL;
}

//----------------------------------------------------------------------------
//Destructor
vtkSlicerPlannerLogic::~vtkSlicerPlannerLogic()
{
}

//-----------------------------------------------------------------------------
const char* vtkSlicerPlannerLogic::DeleteChildrenWarningSettingName()
{
  return "Planner/DeleteChildrenWarning";
}

//----------------------------------------------------------------------------
bool vtkSlicerPlannerLogic::DeleteHierarchyChildren(vtkMRMLNode* node)
{
  vtkMRMLHierarchyNode* hNode = vtkMRMLHierarchyNode::SafeDownCast(node);
  if(!hNode)
  {
    vtkErrorMacro("DeleteHierarchyChildren: Not a hierarchy node.");
    return false;
  }
  if(!this->GetMRMLScene())
  {
    vtkErrorMacro("DeleteHierarchyChildren: No scene defined on this class");
    return false;
  }

  // first off, set up batch processing mode on the scene
  this->GetMRMLScene()->StartState(vtkMRMLScene::BatchProcessState);

  // get all the children nodes
  std::vector< vtkMRMLHierarchyNode*> allChildren;
  hNode->GetAllChildrenNodes(allChildren);

  // and loop over them
  for(unsigned int i = 0; i < allChildren.size(); ++i)
  {
    vtkMRMLNode* associatedNode = allChildren[i]->GetAssociatedNode();
    if(associatedNode)
    {
      this->GetMRMLScene()->RemoveNode(associatedNode);
    }
    this->GetMRMLScene()->RemoveNode(allChildren[i]);
  }
  // end batch processing
  this->GetMRMLScene()->EndState(vtkMRMLScene::BatchProcessState);

  return true;
}

//----------------------------------------------------------------------------
void vtkSlicerPlannerLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
void vtkSlicerPlannerLogic::SetMRMLSceneInternal(vtkMRMLScene* newScene)
{
  Superclass::SetMRMLSceneInternal(newScene);
}

//---------------------------------------------------------------------------
void vtkSlicerPlannerLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//----------------------------------------------------------------------------
//Set logic for Shrink Wrap CLI
void vtkSlicerPlannerLogic::setWrapperLogic(vtkSlicerCLIModuleLogic* logic)
{
  this->wrapperLogic = logic;
}

//----------------------------------------------------------------------------
//Create reference model form current hierarhcy state
vtkMRMLCommandLineModuleNode* vtkSlicerPlannerLogic::createPreOPModels(vtkMRMLModelHierarchyNode* HierarchyNode)
{
  if(this->SkullWrappedPreOP)
  {
    this->GetMRMLScene()->RemoveNode(this->SkullWrappedPreOP);
    this->SkullWrappedPreOP = NULL;
  }

  std::string name;
  name = HierarchyNode->GetName();
  name += " - Merged";
  this->TempMerged = this->mergeModel(HierarchyNode, name);
  this->TempMerged->GetDisplayNode()->SetVisibility(0);
  name = HierarchyNode->GetName();
  name += " - Wrapped";
  this->saveModelHierarchyAsMRMLScene(HierarchyNode, "irrelevant");
  return this->wrapModel(this->TempMerged, name, vtkSlicerPlannerLogic::PreOP);
}

//----------------------------------------------------------------------------
//Get the pre-op ICV
double vtkSlicerPlannerLogic::getPreOPICV()
{
  if(this->SkullWrappedPreOP)
  {
    this->preOPICV = this->computeICV(this->SkullWrappedPreOP);
  }

  return this->preOPICV;
}

//----------------------------------------------------------------------------
//Create wrapped model from current hierarchy
vtkMRMLCommandLineModuleNode*  vtkSlicerPlannerLogic::createCurrentModel(vtkMRMLModelHierarchyNode* HierarchyNode)
{
  if(this->CurrentModel)
  {
    this->GetMRMLScene()->RemoveNode(this->CurrentModel);
    this->CurrentModel = NULL;
  }
  std::string name;
  name = HierarchyNode->GetName();
  name += " - Temp Merge";
  this->TempMerged = this->mergeModel(HierarchyNode, name);
  this->TempMerged->GetDisplayNode()->SetVisibility(0);
  name = HierarchyNode->GetName();
  name += " - Current Wrapped";
  return this->wrapModel(this->TempMerged, name, vtkSlicerPlannerLogic::Current);
}

//----------------------------------------------------------------------------
//Get the current ICV
double vtkSlicerPlannerLogic::getCurrentICV()
{
  if(this->CurrentModel)
  {
    this->currentICV = this->computeICV(this->CurrentModel);
  }
  return this->currentICV;
}

//----------------------------------------------------------------------------
//Create wrapped version of brain model input
vtkMRMLCommandLineModuleNode* vtkSlicerPlannerLogic::createHealthyBrainModel(vtkMRMLModelNode* model)
{
  if(this->HealthyBrain)
  {
    this->GetMRMLScene()->RemoveNode(this->HealthyBrain);
    this->HealthyBrain = NULL;
  }

  std::string name;
  name = model->GetName();
  name += " - Wrapped";
  return wrapModel(model, name, vtkSlicerPlannerLogic::Brain);
}

//----------------------------------------------------------------------------
//Get brain ICV
double vtkSlicerPlannerLogic::getHealthyBrainICV()
{
  if(this->HealthyBrain)
  {
    this->healthyBrainICV = this->computeICV(this->HealthyBrain);
  }
  return this->healthyBrainICV;
}

//----------------------------------------------------------------------------
//Create wrapped version of bone template input
vtkMRMLCommandLineModuleNode* vtkSlicerPlannerLogic::createBoneTemplateModel(vtkMRMLModelNode* model)
{
  if (this->BoneTemplate)
  {
    this->GetMRMLScene()->RemoveNode(this->BoneTemplate);
    this->BoneTemplate = NULL;
  }

  std::string name;
  name = model->GetName();
  name += " - Wrapped";
  return wrapModel(model, name, vtkSlicerPlannerLogic::Template);
}

//----------------------------------------------------------------------------
//Get template ICV
double vtkSlicerPlannerLogic::getTemplateICV()
{
  if (this->BoneTemplate)
  {
    this->templateICV = this->computeICV(this->BoneTemplate);
  }
  return this->templateICV;
}

//----------------------------------------------------------------------------
//Merge hierarchy into a single model
vtkMRMLModelNode* vtkSlicerPlannerLogic::mergeModel(vtkMRMLModelHierarchyNode* HierarchyNode, std::string name)
{

  vtkNew<vtkMRMLModelNode> mergedModel;
  vtkNew<vtkAppendPolyData> filter;
  vtkMRMLScene* scene = this->GetMRMLScene();
  mergedModel->SetScene(scene);
  mergedModel->SetName(name.c_str());
  vtkNew<vtkMRMLModelDisplayNode> dnode;
  vtkNew<vtkMRMLModelStorageNode> snode;
  mergedModel->SetAndObserveDisplayNodeID(dnode->GetID());
  mergedModel->SetAndObserveStorageNodeID(snode->GetID());
  scene->AddNode(dnode.GetPointer());
  scene->AddNode(snode.GetPointer());
  scene->AddNode(mergedModel.GetPointer());

  std::vector<vtkMRMLHierarchyNode*> children;
  std::vector<vtkMRMLHierarchyNode*>::const_iterator it;
  HierarchyNode->GetAllChildrenNodes(children);
  for(it = children.begin(); it != children.end(); ++it)
  {
    vtkMRMLModelNode* childModel =
      vtkMRMLModelNode::SafeDownCast((*it)->GetAssociatedNode());

    if(childModel)
    {
      filter->AddInputData(childModel->GetPolyData());

    }
  }

  filter->Update();
  mergedModel->SetAndObservePolyData(filter->GetOutput());
  mergedModel->SetAndObserveDisplayNodeID(dnode->GetID());

  return mergedModel.GetPointer();
}

//----------------------------------------------------------------------------
//Compute the ICV of a model
double vtkSlicerPlannerLogic::computeICV(vtkMRMLModelNode* model)
{
  vtkNew<vtkMassProperties> areaFilter;
  vtkNew<vtkTriangleFilter> triFilter;
  triFilter->SetInputData(model->GetPolyData());
  triFilter->Update();
  areaFilter->SetInputData(triFilter->GetOutput());
  areaFilter->Update();
  return (areaFilter->GetVolume() / 1000);   //convert to cm^3
}

//----------------------------------------------------------------------------
//Create shrink wrapped version of a model
vtkMRMLCommandLineModuleNode* vtkSlicerPlannerLogic::wrapModel(vtkMRMLModelNode* model, std::string name, int dest)
{
  vtkNew<vtkMRMLModelNode> wrappedModel;
  vtkMRMLScene* scene = this->GetMRMLScene();
  wrappedModel->SetScene(scene);
  wrappedModel->SetName(name.c_str());
  vtkNew<vtkMRMLModelDisplayNode> dnode;
  vtkNew<vtkMRMLModelStorageNode> snode;
  wrappedModel->SetAndObserveDisplayNodeID(dnode->GetID());
  wrappedModel->SetAndObserveStorageNodeID(snode->GetID());
  scene->AddNode(dnode.GetPointer());
  scene->AddNode(snode.GetPointer());
  scene->AddNode(wrappedModel.GetPointer());

  switch(dest)
  {
  case vtkSlicerPlannerLogic::Current:
    this->CurrentModel = wrappedModel.GetPointer();
    break;
  case vtkSlicerPlannerLogic::PreOP:
    this->SkullWrappedPreOP = wrappedModel.GetPointer();
    break;
  case vtkSlicerPlannerLogic::Brain:
    this->HealthyBrain = wrappedModel.GetPointer();
    break;
  case vtkSlicerPlannerLogic::Template:
    this->BoneTemplate = wrappedModel.GetPointer();

  }

  //CLI setup
  this->wrapperLogic->SetMRMLScene(this->GetMRMLScene());
  vtkMRMLCommandLineModuleNode* cmdNode = this->wrapperLogic->CreateNodeInScene();
  cmdNode->SetParameterAsString("inputModel", model->GetID());
  cmdNode->SetParameterAsString("outputModel", wrappedModel->GetID());
  cmdNode->SetParameterAsString("PhiRes", "20");
  cmdNode->SetParameterAsString("ThetaRes", "20");
  this->wrapperLogic->Apply(cmdNode, true);
  return cmdNode;
}

//----------------------------------------------------------------------------
//Finish up wrapper CLI
void vtkSlicerPlannerLogic::finishWrap(vtkMRMLCommandLineModuleNode* cmdNode)
{
  vtkMRMLModelNode* node = vtkMRMLModelNode::SafeDownCast(this->GetMRMLScene()->GetNodeByID(cmdNode->GetParameterAsString("outputModel")));
  node->GetDisplayNode()->SetVisibility(0);
  node->GetDisplayNode()->SetActiveScalarName("Normals");
  this->GetMRMLScene()->RemoveNode(cmdNode);
  node->HideFromEditorsOn();
  node->SetAttribute("PlannerRole", "WrappedModel");

  if(this->TempMerged)
  {
    this->GetMRMLScene()->RemoveNode(this->TempMerged);
    this->TempMerged = NULL;
  }

  if(this->TempWrapped)
  {
    this->GetMRMLScene()->RemoveNode(this->TempWrapped);
    this->TempWrapped = NULL;
  }
}

//----------------------------------------------------------------------------
//Fill table node with metrics
void vtkSlicerPlannerLogic::fillMetricsTable(vtkMRMLModelHierarchyNode* HierarchyNode, vtkMRMLTableNode* modelMetricsTable)
{
  double preOpVolume;
  double currentVolume;
  double brainVolume;
  double templateVolume;
  if(HierarchyNode)
  {
    preOpVolume = this->getPreOPICV();
    brainVolume = this->getHealthyBrainICV();
    currentVolume = this->getCurrentICV();
    templateVolume = this->getTemplateICV();

    modelMetricsTable->RemoveAllColumns();
    std::string modelTableName = "Model Metrics - ";
    modelTableName += HierarchyNode->GetName();
    modelMetricsTable->SetName(modelTableName.c_str());

    modelMetricsTable->AddColumn();
    vtkAbstractArray* col1 = modelMetricsTable->AddColumn();
    col1->SetName("Healthy Brain");
    vtkAbstractArray* col2 = modelMetricsTable->AddColumn();
    col2->SetName("Bone Template");
    vtkAbstractArray* col3 = modelMetricsTable->AddColumn();
    col3->SetName("Pre-op");
    vtkAbstractArray* col4 = modelMetricsTable->AddColumn();
    col4->SetName("Current");
    modelMetricsTable->SetUseColumnNameAsColumnHeader(true);
    modelMetricsTable->SetUseFirstColumnAsRowHeader(true);
    modelMetricsTable->SetLocked(true);

    modelMetricsTable->AddEmptyRow();
    modelMetricsTable->SetCellText(0, 0, "ICV\n cm^3");
    
    std::stringstream brainVolumeSstr;
    brainVolumeSstr << brainVolume;
    const std::string& brainVolumeString = brainVolumeSstr.str();
    modelMetricsTable->SetCellText(0, 1, brainVolumeString.c_str());

    std::stringstream templateVolumeSstr;
    templateVolumeSstr << templateVolume;
    const std::string& templateVolumeString = templateVolumeSstr.str();
    modelMetricsTable->SetCellText(0, 2, templateVolumeString.c_str());
    
    std::stringstream preOpVolumeSstr;
    preOpVolumeSstr << preOpVolume;
    const std::string& preOpVolumeString = preOpVolumeSstr.str();
    modelMetricsTable->SetCellText(0, 3, preOpVolumeString.c_str());
    
    std::stringstream currentVolumeSstr;
    currentVolumeSstr << currentVolume;
    const std::string& currentVolumeString = currentVolumeSstr.str();
    modelMetricsTable->SetCellText(0, 4, currentVolumeString.c_str());
    
  }
}

//----------------------------------------------------------------------------
//Initiaize bending
void vtkSlicerPlannerLogic::initializeBend(vtkPoints* inputFiducials, vtkMRMLModelNode* model)
{
  this->Fiducials = inputFiducials;
  this->ModelToBend = model;

  vtkNew<vtkTriangleFilter> triangulate;
  vtkNew<vtkCleanPolyData> clean;
  vtkNew<vtkPolyDataNormals> normals;
  normals->SetComputePointNormals(1);
  normals->SetComputeCellNormals(1);
  normals->SetAutoOrientNormals(1);
  normals->SetInputData(this->ModelToBend->GetPolyData());
  normals->Update();
  clean->SetInputData(normals->GetOutput());
  clean->Update();
  this->BendingPolyData = clean->GetOutput();

  this->cellLocator = vtkSmartPointer<vtkCellLocator>::New();
  this->cellLocator->SetDataSet(this->BendingPolyData);
  this->cellLocator->BuildLocator();

  this->generateSourcePoints();
  this->bendInitialized =  true;
}

//----------------------------------------------------------------------------
//CReate bend transform based on points and bend magnitude
vtkSmartPointer<vtkThinPlateSplineTransform> vtkSlicerPlannerLogic::getBendTransform(double magnitude)
{
  vtkSmartPointer<vtkThinPlateSplineTransform> transform = vtkSmartPointer<vtkThinPlateSplineTransform>::New();
  if(this->bendInitialized)
  {

    this->TargetPoints = vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < this->SourcePointsDense->GetNumberOfPoints(); i++)
    {
      double p[3];
      this->SourcePointsDense->GetPoint(i, p);
      vtkVector3d point = (vtkVector3d)p;
      vtkVector3d bent = point;
      if(this->bendMode == Double)
      {
        bent = this->bendPoint(point, magnitude);
      }
      if(this->bendMode == Single)
      {
        if(this->bendSide == A)
        {
          if(this->BendingPlane->EvaluateFunction(point.GetData())*this->BendingPlane->EvaluateFunction(this->SourcePoints->GetPoint(0)) > 0)
          {
            bent = this->bendPoint(point, magnitude);
          }
        }
        if(this->bendSide == B)
        {
          if(this->BendingPlane->EvaluateFunction(point.GetData())*this->BendingPlane->EvaluateFunction(this->SourcePoints->GetPoint(1)) > 0)
          {
            bent = this->bendPoint(point, magnitude);
          }
        }

      }
      this->TargetPoints->InsertPoint(i, bent.GetData());
    }

    transform->SetSigma(.0001);
    transform->SetBasisToR();
    transform->SetSourceLandmarks(this->SourcePointsDense);
    transform->SetTargetLandmarks(this->TargetPoints);
    transform->Update();
  }
  return transform;
}

//----------------------------------------------------------------------------
//Clear all bending data
void vtkSlicerPlannerLogic::clearBendingData()
{
  this->SourcePoints = NULL;
  this->SourcePointsDense = NULL;
  this->TargetPoints = NULL;
  this->Fiducials = NULL;
  this->ModelToBend = NULL;
  this->cellLocator = NULL;
  this->BendingPlane = NULL;
  this->BendingPlaneLocator = NULL;
  this->bendInitialized = false;
}

//----------------------------------------------------------------------------
//Create source points based on fiducials
void vtkSlicerPlannerLogic::generateSourcePoints()
{
  this->SourcePoints = vtkSmartPointer<vtkPoints>::New();
  bool bendingAxis = true;

  //A and B are on the bending line.  C and D are on the bending axis
  vtkVector3d A;
  vtkVector3d B;
  vtkVector3d C;
  vtkVector3d D;
  vtkVector3d CD;
  vtkVector3d AB;

  
  double firstPoint[3];
  double secondPoint[3];

  this->Fiducials->GetPoint(0, firstPoint);
  this->Fiducials->GetPoint(1, secondPoint);

  //we are receive the bending axis as input, must derive bending line
  if (bendingAxis)
  {
      C = (vtkVector3d)firstPoint;
      D = (vtkVector3d)secondPoint;
      CD = D - C;
      vtkVector3d CDMid = C + 0.5*CD;
      vtkVector3d normal = this->getNormalAtPoint(CDMid, this->cellLocator, this->BendingPolyData);
      vtkVector3d bendLine = normal.Cross(CD);
      A = CDMid + bendLine;
      B = CDMid - bendLine;
      AB = B - A;
  }
  else  //we are receiving the bending line as input, must derive bending axis
  {
      A = (vtkVector3d)firstPoint;
      B = (vtkVector3d)secondPoint;
      AB = B - A;
      vtkVector3d ABMid = A + 0.5*AB;
      vtkVector3d normal = this->getNormalAtPoint(ABMid, this->cellLocator, this->BendingPolyData);
      vtkVector3d bendAxis = normal.Cross(AB);
      C = ABMid + bendAxis;
      D = ABMid - bendAxis;
      CD = D - C;
  }
    
  AB.Normalize();
  CD.Normalize();

  vtkSmartPointer<vtkPlane> fixedPlane = this->createPlane(C, D, A, B);
  this->BendingPlane = fixedPlane;
  this->createBendingLocator();

  
  A = this->projectToModel(A);
  B = this->projectToModel(B);

  this->SourcePoints->InsertPoint(0, A.GetData());
  this->SourcePoints->InsertPoint(1, B.GetData());
  this->SourcePoints->InsertPoint(2, C.GetData());
  this->SourcePoints->InsertPoint(3, D.GetData());

  //Compute Next point as the vector defining the bend axis
  vtkVector3d E = A + 0.5 * (B - A);
  vtkVector3d CE = E - C;
  CD = D - C;
  vtkVector3d CF = CE.Dot(CD.Normalized()) * CD.Normalized();

  //Midpoint projected onto te line between the fixed points - Pivot point
  vtkVector3d F = C + CF;
  F = projectToModel(F);
  vtkVector3d FE = E - F;
  vtkVector3d FB = B - F;
  vtkVector3d axis = FE.Cross(FB);
  axis.Normalize();

  //Store beding axis in source points
  this->SourcePoints->InsertPoint(4, axis.GetData());

  //Agressively downsample to create source points
  vtkNew<vtkCleanPolyData> clean;
  vtkNew<vtkVertexGlyphFilter> verts;
  verts->SetInputData(this->BendingPolyData);
  verts->Update();
  clean->SetInputData(verts->GetOutput());
  clean->SetTolerance(0.07);
  clean->Update();
  this->SourcePointsDense = clean->GetOutput()->GetPoints();
}

//----------------------------------------------------------------------------
//Project a 3D point onto the closest point on the bending model
vtkVector3d vtkSlicerPlannerLogic::projectToModel(vtkVector3d point)
{
  //build locator when model is loaded

  return this->projectToModel(point, this->cellLocator);
}

//----------------------------------------------------------------------------
//Project a 3D point onto the closest point on the bending model, constrained by a plane
vtkVector3d vtkSlicerPlannerLogic::projectToModel(vtkVector3d point, vtkPlane* plane)
{
  vtkNew<vtkCutter> cutter;
  cutter->SetCutFunction(plane);
  cutter->SetInputData(this->BendingPolyData);
  vtkSmartPointer<vtkPolyData> cut;
  cutter->Update();
  cut = cutter->GetOutput();
  return this->projectToModel(point, cut);
}

//----------------------------------------------------------------------------
//Project a 3D point onto the closest point on the specified model
vtkVector3d vtkSlicerPlannerLogic::projectToModel(vtkVector3d point, vtkPolyData* model)
{
  vtkNew<vtkCellLocator> locator;
  vtkNew<vtkTriangleFilter> triangulate;
  triangulate->SetInputData(model);
  triangulate->Update();
  locator->SetDataSet(triangulate->GetOutput());
  locator->BuildLocator();
  return this->projectToModel(point, locator.GetPointer());
}

//----------------------------------------------------------------------------
//Project a 3D point onto the closest point on the model as defined by the provided cell locator
vtkVector3d vtkSlicerPlannerLogic::projectToModel(vtkVector3d point, vtkCellLocator* locator)
{
  double closestPoint[3];//the coordinates of the closest point will be returned here
  double closestPointDist2; //the squared distance to the closest point will be returned here
  vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
  int subId; //this is rarely used (in triangle strips only, I believe)
  locator->FindClosestPoint(point.GetData(), closestPoint, cellId, subId, closestPointDist2);

  vtkVector3d projection;
  projection.Set(closestPoint[0], closestPoint[1], closestPoint[2]);
  return projection;
}

//----------------------------------------------------------------------------
//Create Plane from two points in plane and two points on normal vector
vtkSmartPointer<vtkPlane> vtkSlicerPlannerLogic::createPlane(vtkVector3d A, vtkVector3d B, vtkVector3d C, vtkVector3d D)
{
  //A and B are in the plane
  //C and D are perp to plane

  vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
  vtkVector3d AB = B - A;
  vtkVector3d E = A + 0.5 * AB;
  vtkVector3d CD = D - C;
  plane->SetOrigin(E.GetData());
  plane->SetNormal(CD.GetData());
  return plane;
}

//----------------------------------------------------------------------------
//bend point using vector computed from axis
vtkVector3d vtkSlicerPlannerLogic::bendPoint(vtkVector3d point, double magnitude)
{
  double ax[3];
  this->SourcePoints->GetPoint(4, ax);
  vtkVector3d axis = (vtkVector3d)ax;
  vtkVector3d F = projectToModel(point, this->BendingPlaneLocator);
  vtkVector3d AF = F - point;
  vtkVector3d BendingVector;
  if(this->BendingPlane->EvaluateFunction(point.GetData()) < 0)
  {
    BendingVector = AF.Cross(axis);
  }
  else
  {
    BendingVector = axis.Cross(AF);
  }
  vtkVector3d point2 = point + ((magnitude * AF.Norm()) * BendingVector.Normalized());
  vtkVector3d A2F = F - point2;

  //correction factor
  point2 = point2 + (A2F.Norm() - AF.Norm()) * A2F.Normalized();

  return point2;
}
//----------------------------------------------------------------------------
//Create a point locator constrained to the bending axis
void vtkSlicerPlannerLogic::createBendingLocator()
{
  this->BendingPlaneLocator = vtkSmartPointer<vtkCellLocator>::New();

  vtkNew<vtkCutter> cutter;
  cutter->SetCutFunction(this->BendingPlane);
  cutter->SetInputData(this->BendingPolyData);
  vtkSmartPointer<vtkPolyData> cut;
  cutter->Update();
  cut = cutter->GetOutput();

  vtkNew<vtkTriangleFilter> triangulate;
  triangulate->SetInputData(cut);
  triangulate->Update();
  this->BendingPlaneLocator->SetDataSet(triangulate->GetOutput());
  this->BendingPlaneLocator->BuildLocator();
}

//----------------------------------------------------------------------------
//Remove models and clear data
void vtkSlicerPlannerLogic::clearModelsAndData()
{
  this->clearBendingData();
  if (this->SkullWrappedPreOP)
  {
    this->GetMRMLScene()->RemoveNode(this->SkullWrappedPreOP);
    this->SkullWrappedPreOP = NULL;
  }
  if (this->HealthyBrain)
  {
    this->GetMRMLScene()->RemoveNode(this->HealthyBrain);
    this->HealthyBrain = NULL;
  }
  if (this->CurrentModel)
  {
    this->GetMRMLScene()->RemoveNode(this->CurrentModel);
    this->CurrentModel = NULL;
  }
  if (this->BoneTemplate)
  {
    this->GetMRMLScene()->RemoveNode(this->BoneTemplate);
    this->BoneTemplate = NULL;
  }
  if (this->TempMerged)
  {
    this->GetMRMLScene()->RemoveNode(this->TempMerged);
    this->TempMerged = NULL;
  }

  if (this->TempWrapped)
  {
    this->GetMRMLScene()->RemoveNode(this->TempWrapped);
    this->TempWrapped = NULL;
  }

  this->preOPICV = 0;
  this->healthyBrainICV = 0;
  this->currentICV = 0;
  this->templateICV = 0;

}

vtkVector3d vtkSlicerPlannerLogic::getNormalAtPoint(vtkVector3d point, vtkCellLocator* locator, vtkPolyData* model)
{
    double closestPoint[3];//the coordinates of the closest point will be returned here
    double closestPointDist2; //the squared distance to the closest point will be returned here
    vtkIdType cellId; //the cell id of the cell containing the closest point will be returned here
    int subId; //this is rarely used (in triangle strips only, I believe)
    locator->FindClosestPoint(point.GetData(), closestPoint, cellId, subId, closestPointDist2);

    vtkVector3d normal;
    double n[3];
    model->GetCellData()->GetNormals()->GetTuple(cellId, n);
    normal.SetX(n[0]);
    normal.SetY(n[1]);
    normal.SetZ(n[2]);

    return normal;
}

void vtkSlicerPlannerLogic::AddCompleteModelHierarchyToMiniScene(vtkMRMLScene *miniscene, vtkMRMLModelHierarchyNode *mhnd)
{
    typedef std::map<std::string, std::string> IDMapType;
    IDMapType IDMap;
    std::string PlanDirectory = "D:/Plans/TestPlan";
    if (mhnd)
    {
        // construct a list that includes this node and all its children
        std::vector<vtkMRMLHierarchyNode*> hnodes;
        mhnd->GetAllChildrenNodes(hnodes);
        hnodes.insert(hnodes.begin(), mhnd);  // add the current node to the front of the vector

                                              // copy the entire hierarchy into the miniscene, we assume the nodes are ordered such that parents appear before children
        for (std::vector<vtkMRMLHierarchyNode*>::iterator it = hnodes.begin(); it != hnodes.end(); ++it)
        {
            vtkMRMLNode *tnd = *it;
            vtkMRMLModelHierarchyNode *tmhnd = vtkMRMLModelHierarchyNode::SafeDownCast(tnd);

            if (!tmhnd)
            {
                std::cerr << "Child is not a model hierarchy node." << std::endl;
                continue;
            }

            // model hierarchy nodes need to get put in a scene
            vtkMRMLNode *cp = miniscene->CopyNode(tnd);
            vtkMRMLModelHierarchyNode *mhcp = vtkMRMLModelHierarchyNode::SafeDownCast(cp);

            // wire the parent relationship (again, we assume the parents appeared in the list before the children)
            vtkMRMLNode *p = tmhnd->GetParentNode();
            if (p)
            {
                // find parent in the sceneToMiniSceneMap
                IDMapType::iterator mit = IDMap.find(p->GetID());
                if (mit != IDMap.end())
                {
                    mhcp->SetParentNodeID((*mit).second.c_str());
                }
            }

            // keep track of what scene node corresponds to what miniscene node
            IDMap[tnd->GetID()] = cp->GetID();

            // also add any display node
            vtkMRMLDisplayNode *dnd = tmhnd->GetDisplayNode();
            if (dnd)
            {
                vtkMRMLNode *dcp = miniscene->CopyNode(dnd);

                vtkMRMLDisplayNode *d = vtkMRMLDisplayNode::SafeDownCast(dcp);

                mhcp->SetAndObserveDisplayNodeID(d->GetID());
            }

            // add the actual model node
            vtkMRMLDisplayableNode* tmnd = tmhnd->GetDisplayableNode();
            if (tmnd)
            {
                vtkMRMLModelNode *mnd = vtkMRMLModelNode::SafeDownCast(tmnd);
                if (mnd)
                {
                    vtkMRMLNode *mcp = miniscene->CopyNode(mnd);
                    vtkMRMLModelNode *tmcp = vtkMRMLModelNode::SafeDownCast(mcp);

                    mhcp->SetAssociatedNodeID(tmcp->GetID());

                    // add the display nodes for the model to the miniscene
                    int ndn = mnd->GetNumberOfDisplayNodes();
                    for (int i = 0; i<ndn; i++)
                    {
                        vtkMRMLDisplayNode *mdnd = mnd->GetNthDisplayNode(i);
                        if (mdnd && tmcp)
                        {
                            vtkMRMLNode *mdcp = miniscene->CopyNode(mdnd);

                            vtkMRMLDisplayNode *d = vtkMRMLDisplayNode::SafeDownCast(mdcp);

                            tmcp->SetAndObserveDisplayNodeID(d->GetID());
                        }
                    }

                    // add the storage node for the model to the miniscene
                    vtkMRMLStorageNode *msnd = mnd->GetStorageNode();
                    if (msnd)
                    {
                        vtkMRMLNode *mscp = miniscene->CopyNode(msnd);

                        vtkMRMLModelStorageNode *s = vtkMRMLModelStorageNode::SafeDownCast(mscp);
                        std::string fname = PlanDirectory + "/" + tmcp->GetName() + "_" + tmcp->GetID() + ".vtk";
                        s->SetFileName(fname.c_str());
                        if (tmcp)
                        {
                            tmcp->SetAndObserveStorageNodeID(s->GetID());
                        }
                    }

                    // keep track of the what scene node corresponds to what miniscene node
                    IDMap[mnd->GetID()] = mcp->GetID();
                }
            }
        }
    }
}

void vtkSlicerPlannerLogic::saveModelHierarchyAsMRMLScene(vtkMRMLModelHierarchyNode* HierarchyNode, std::string PlanDirectory)
{
    //Create output scene
    vtkNew<vtkMRMLScene> outputScene;
    std::string planDirect = "D:/Plans/TestPlan";
    std::string mrmLFilename = planDirect + "/" + "testScene.mrb";
    outputScene->SetRootDirectory(vtksys::SystemTools::GetParentDirectory(mrmLFilename.c_str()).c_str());

    this->AddCompleteModelHierarchyToMiniScene(outputScene.GetPointer(), HierarchyNode);
    this->writeToMRB(outputScene, mrmLFilename);
    outputScene->Clear(1);
    
}

//----------------------------------------------------------------------------
bool vtkSlicerPlannerLogic::writeToMRB(vtkMRMLScene * scene, std::string filename)
{
    //
    // make a temp directory to save the scene into - this will
    // be a uniquely named directory that contains a directory
    // named based on the user's selection.
    //

    QFileInfo fileInfo(filename.c_str());
    QString basePath = fileInfo.absolutePath();
    if (!QFileInfo(basePath).isWritable())
    {
        qWarning() << "Failed to save" << fileInfo.absoluteFilePath() << ":"
            << "Path" << basePath << "is not writable";
        return false;
    }

    // TODO: switch to QTemporaryDir in Qt5.
    // For now, create a named directory and use Qt calls to remove it
    QString tempDir = qSlicerCoreApplication::application()->temporaryPath();
    QFileInfo pack(QDir(tempDir), //QDir::tempPath(),
        QString("__BundleSaveTemp-") +
        QDateTime::currentDateTime().toString("yyyy-MM-dd_hh+mm+ss.zzz"));
    qDebug() << "packing to " << pack.absoluteFilePath();

    // make a subdirectory with the name the user has chosen
    QFileInfo bundle = QFileInfo(QDir(pack.absoluteFilePath()),
        fileInfo.completeBaseName());
    QString bundlePath = bundle.absoluteFilePath();
    if (bundle.exists())
    {
        if (!ctk::removeDirRecursively(bundlePath))
        {
            QMessageBox::critical(0, "Save Scene as MRB", "Could not remove temp directory");
            return false;
        }
    }

    if (!QDir().mkpath(bundlePath))
    {
        QMessageBox::critical(0, "Save scene as MRB", "Could not make temp directory");
        return false;
    }    

    //
    // Now save the scene into the bundle directory and then make a zip (mrb) file
    // in the user's selected file location
    //
    vtkSlicerApplicationLogic* applicationLogic =
        qSlicerCoreApplication::application()->applicationLogic();
    bool retval =
        this->SaveSceneToSlicerDataBundleDirectory(bundlePath.toLatin1(), scene);
    if (!retval)
    {
        QMessageBox::critical(0, "Save scene as MRB", "Failed to create bundle");
        return false;
    }

    qDebug() << "zipping to " << fileInfo.absoluteFilePath();
    if (!applicationLogic->Zip(fileInfo.absoluteFilePath().toLatin1(),
        bundlePath.toLatin1()))
    {
        QMessageBox::critical(0, "Save scene as MRB", "Could not compress bundle");
        return false;
    }

    //
    // Now clean up the temp directory
    //
    if (!ctk::removeDirRecursively(bundlePath))
    {
        QMessageBox::critical(0, "Save scene as MRB", "Could not remove temp directory");
        return false;
    }

    // Mark the storable nodes as modified since read, since that flag was reset
    // when the files were written out. If there was newly generated data in the
    // scene that only got saved to the MRB bundle directory, it would be marked
    // as unmodified since read when saving as a MRML file + data. This will not
    // disrupt multiple MRB saves.
    scene->SetStorableNodesModifiedSinceRead();

    qDebug() << "saved " << fileInfo.absoluteFilePath();
    return true;
}

//----------------------------------------------------------------------------
bool vtkSlicerPlannerLogic::SaveSceneToSlicerDataBundleDirectory(const char *sdbDir, vtkMRMLScene *scene)
{

    //
    // first, confirm the arguments are valid and create directories if needed
    // then, save all paths and filenames in the current scene
    //  and replace them with paths to the sdbDir
    // then create a scene view of the contents of the data bundle
    // then save the scene
    // -- replace the original paths
    // -- remove the scene view
    //
    // at the end, the scene should be restored to its original state
    // except that some storables will have default storage nodes
    //

    
    if (!sdbDir)
    {
        vtkErrorMacro("no directory given!");
        return false;
    }

    // if the path to the directory is not absolute, return
    if (!vtksys::SystemTools::FileIsFullPath(sdbDir))
    {
        vtkErrorMacro("given directory is not a full path: " << sdbDir);
        return false;
    }
    // is it a directory?
    if (!vtksys::SystemTools::FileIsDirectory(sdbDir))
    {
        vtkErrorMacro("given directory name is not actually a directory, try again!" << sdbDir);
        return false;
    }
    std::string rootDir = std::string(sdbDir);
    vtkDebugMacro("Using root dir of " << rootDir);
    // need the components to build file names
    std::vector<std::string> rootPathComponents;
    vtksys::SystemTools::SplitPath(rootDir.c_str(), rootPathComponents);

    // remove the directory if it does exist
    if (vtksys::SystemTools::FileExists(rootDir.c_str(), false))
    {
        if (!vtksys::SystemTools::RemoveADirectory(rootDir.c_str()))
        {
            vtkErrorMacro("Error removing SDB scene directory " << rootDir.c_str() << ", cannot make a fresh archive.");
            return false;
        }
    }
    // create the SDB directory
    if (!vtksys::SystemTools::FileExists(rootDir.c_str(), false))
    {
        if (!vtksys::SystemTools::MakeDirectory(rootDir.c_str()))
        {
            vtkErrorMacro("Unable to make temporary directory " << rootDir);
            return false;
        }
    }

    //
    // now, replace paths with data bundle paths, saving the original values
    //

    // the root directory
    std::string origURL(scene->GetURL());
    std::string origRootDirectory(scene->GetRootDirectory());

    // the new url of the mrml scene
    std::string urlStr = vtksys::SystemTools::GetFilenameWithoutExtension(rootDir.c_str()) + std::string(".mrml");
    rootPathComponents.push_back(urlStr);
    urlStr = vtksys::SystemTools::JoinPath(rootPathComponents);
    rootPathComponents.pop_back();
    vtkDebugMacro("set new scene url to " << scene->GetURL());

    // the new data directory
    std::vector<std::string> pathComponents;
    vtksys::SystemTools::SplitPath(rootDir.c_str(), pathComponents);
    pathComponents.push_back("Data");
    std::string dataDir = vtksys::SystemTools::JoinPath(pathComponents);
    vtkDebugMacro("using data dir of " << dataDir);

    // create the data dir
    if (!vtksys::SystemTools::FileExists(dataDir.c_str()))
    {
        if (!vtksys::SystemTools::MakeDirectory(dataDir.c_str()))
        {
            vtkErrorMacro("Unable to make data directory " << dataDir);
            return false;
        }
    }

    //
    // start changing the scene - don't return from below here
    // until scene has been restored to original state
    //
    scene->SetRootDirectory(rootDir.c_str());
    scene->SetURL(urlStr.c_str());

    // change all storage nodes and file names to be unique in the new directory
    // write the new data as we go; save old values
    this->OriginalStorageNodeDirs.clear();
    this->OriginalStorageNodeFileNames.clear();

    std::map<std::string, vtkMRMLNode *> storableNodes;

    int numNodes = scene->GetNumberOfNodes();
    for (int i = 0; i < numNodes; ++i)
    {
        vtkMRMLNode *mrmlNode = scene->GetNthNode(i);
        if (!mrmlNode)
        {
            vtkErrorMacro("unable to get " << i << "th node from scene with " << numNodes << " nodes");
            continue;
        }
        if (mrmlNode->IsA("vtkMRMLStorableNode"))
        {
            // get all storable nodes in the main scene
            // and store them in the map by ID to avoid duplicates for the scene views
            vtkMRMLStorableNode *storableNode = vtkMRMLStorableNode::SafeDownCast(mrmlNode);

            this->SaveStorableNodeToSlicerDataBundleDirectory(storableNode, dataDir);

            storableNodes[std::string(storableNode->GetID())] = storableNode;
        }
        if (mrmlNode->IsA("vtkMRMLSceneViewNode"))
        {
            // get all additional storable nodes for all scene views
            vtkMRMLSceneViewNode *sceneViewNode = vtkMRMLSceneViewNode::SafeDownCast(mrmlNode);
            sceneViewNode->SetSceneViewRootDir(scene->GetRootDirectory());

            std::vector<vtkMRMLNode *> snodes;
            sceneViewNode->GetNodesByClass("vtkMRMLStorableNode", snodes);
            std::vector<vtkMRMLNode *>::iterator sit;
            for (sit = snodes.begin(); sit != snodes.end(); sit++)
            {
                vtkMRMLStorableNode* storableNode = vtkMRMLStorableNode::SafeDownCast(*sit);
                if (storableNodes.find(std::string(storableNode->GetID())) == storableNodes.end())
                {
                    // save only new storable nodes
                    storableNode->SetAddToScene(1);
                    storableNode->UpdateScene(scene);
                    this->SaveStorableNodeToSlicerDataBundleDirectory(storableNode, dataDir);

                    storableNodes[std::string(storableNode->GetID())] = storableNode;
                    storableNode->SetAddToScene(0);
                }
                else
                {
                    // just do the path save/update since the paths are indexed by the node, not id
                    vtkMRMLStorageNode *storageNode = storableNode->GetStorageNode();
                    if (storageNode)
                    {
                        std::string fileName(storageNode->GetFileName());
                        this->OriginalStorageNodeFileNames[storageNode].push_back(fileName);
                        for (int i = 0; i < storageNode->GetNumberOfFileNames(); ++i)
                        {
                            this->OriginalStorageNodeFileNames[storageNode].push_back(storageNode->GetNthFileName(i));
                        }
                    }
                }
            }
        }
    }
    //
    // create a scene view, using the snapshot passed in if any
    //
    vtkNew<vtkMRMLSceneViewNode> newSceneViewNode;
    newSceneViewNode->SetScene(scene);
    newSceneViewNode->SetName(scene->GetUniqueNameByString("Slicer Data Bundle Scene View"));
    newSceneViewNode->SetSceneViewDescription("Scene at MRML file save point");
    // save the scene view
    newSceneViewNode->StoreScene();
    scene->AddNode(newSceneViewNode.GetPointer());

    vtkSmartPointer<vtkMRMLStorageNode> newSceneViewStorageNode;
    

    // write the scene to disk, changes paths to relative
    vtkDebugMacro("calling commit on the scene, to url " << scene->GetURL());
    scene->Commit();

    //
    // Now, restore the state of the scene
    //

    scene->SetURL(origURL.c_str());
    scene->SetRootDirectory(origRootDirectory.c_str());

    // clean up scene views
    scene->RemoveNode(newSceneViewNode.GetPointer());
    if (newSceneViewStorageNode)
    {
        scene->RemoveNode(newSceneViewStorageNode);
    }

    // reset the storage paths
    numNodes = scene->GetNumberOfNodes();
    for (int i = 0; i < numNodes; ++i)
    {
        vtkMRMLNode *mrmlNode = scene->GetNthNode(i);
        if (!mrmlNode)
        {
            vtkErrorMacro("unable to get " << i << "th node from scene with " << numNodes << " nodes");
            continue;
        }
        if (mrmlNode->IsA("vtkMRMLSceneViewNode"))
        {
            vtkMRMLSceneViewNode *sceneViewNode = vtkMRMLSceneViewNode::SafeDownCast(mrmlNode);
            sceneViewNode->GetScene()->SetURL(origURL.c_str());
            sceneViewNode->SetSceneViewRootDir(origRootDirectory.c_str());

            // get all additional storable nodes for all scene views
            std::vector<vtkMRMLNode *> snodes;
            sceneViewNode->GetNodesByClass("vtkMRMLStorableNode", snodes);
            std::vector<vtkMRMLNode *>::iterator sit;
            for (sit = snodes.begin(); sit != snodes.end(); sit++)
            {
                vtkMRMLStorableNode* storableNode = vtkMRMLStorableNode::SafeDownCast(*sit);
                vtkMRMLStorageNode *storageNode = storableNode->GetStorageNode();

                if (storageNode && this->OriginalStorageNodeFileNames.find(storageNode) != this->OriginalStorageNodeFileNames.end())
                {
                    storageNode->SetFileName(this->OriginalStorageNodeFileNames[storageNode][0].c_str());
                    if (this->OriginalStorageNodeFileNames[storageNode].size() > 1)
                    {
                        // set the file list
                        storageNode->ResetFileNameList();
                        for (unsigned int fileNumber = 0;
                            fileNumber < this->OriginalStorageNodeFileNames[storageNode].size();
                            ++fileNumber)
                        {
                            // the fileName is also in the file list, but AddFileName does
                            // check for duplicates
                            storageNode->AddFileName(this->OriginalStorageNodeFileNames[storageNode][fileNumber].c_str());
                        }
                    }
                }
                if (storageNode && this->OriginalStorageNodeDirs.find(storageNode) != this->OriginalStorageNodeDirs.end())
                {
                    storageNode->SetDataDirectory(this->OriginalStorageNodeDirs[storageNode].c_str());
                }
            }
        }
        if (mrmlNode->IsA("vtkMRMLStorableNode"))
        {
            vtkMRMLStorableNode *storableNode = vtkMRMLStorableNode::SafeDownCast(mrmlNode);
            vtkMRMLStorageNode *storageNode = storableNode->GetStorageNode();
            if (storageNode && this->OriginalStorageNodeFileNames.find(storageNode) != this->OriginalStorageNodeFileNames.end())
            {
                // std::cout << "Resetting filename on storage node " << storageNode->GetID() << " from " << storageNode->GetFileName() << " back to " << this->OriginalStorageNodeFileNames[storageNode][0].c_str() << "\n\tmodified since read = " << storableNode->GetModifiedSinceRead() << std::endl;
                storageNode->SetFileName(this->OriginalStorageNodeFileNames[storageNode][0].c_str());
                if (this->OriginalStorageNodeFileNames[storageNode].size() > 1)
                {
                    // set the file list
                    for (unsigned int fileNumber = 0;
                        fileNumber < this->OriginalStorageNodeFileNames[storageNode].size();
                        ++fileNumber)
                    {
                        // the fileName is also in the file list
                        storageNode->AddFileName(this->OriginalStorageNodeFileNames[storageNode][fileNumber].c_str());
                    }
                }
            }
            if (storageNode && this->OriginalStorageNodeDirs.find(storageNode) != this->OriginalStorageNodeDirs.end())
            {
                storageNode->SetDataDirectory(this->OriginalStorageNodeDirs[storageNode].c_str());
            }
        }

    }

    scene->SetURL(origURL.c_str());
    scene->SetRootDirectory(origRootDirectory.c_str());

    return true;
}


//-------------------------------------------------------------------------------------------------
void vtkSlicerPlannerLogic::SaveStorableNodeToSlicerDataBundleDirectory(vtkMRMLStorableNode *storableNode,
    std::string &dataDir)
{
    if (!storableNode || !storableNode->GetSaveWithScene())
    {
        return;
    }
    // adjust the file paths for storable nodes
    vtkMRMLStorageNode *storageNode = storableNode->GetStorageNode();
    if (!storageNode)
    {
        vtkDebugMacro("creating a new storage node for " << storableNode->GetID());
        storableNode->AddDefaultStorageNode();
        storageNode = storableNode->GetStorageNode();
        if (!storageNode)
        {
            // no need for storage node to store this node
            return;
        }
    }

    // save the old values for the storage nodes
    // - this->OriginalStorageNodeFileNames has the old filenames (absolute paths)
    // - this->OriginalStorageNodeDirs has old paths
    // std::cout << "SaveStorableNodeToSlicerDataBundleDirectory: saving old storage node file name of " << storageNode->GetFileName() << "\n\tmodified since read = " << storableNode->GetModifiedSinceRead() << std::endl;

    std::string fileName(storageNode->GetFileName() ? storageNode->GetFileName() : "");
    this->OriginalStorageNodeFileNames[storageNode].push_back(fileName);
    for (int i = 0; i < storageNode->GetNumberOfFileNames(); ++i)
    {
        this->OriginalStorageNodeFileNames[storageNode].push_back(storageNode->GetNthFileName(i) ? storageNode->GetNthFileName(i) : "");
    }

    if (fileName.empty())
    {
        // Default storage node usually has empty file name (if Save dialog is not opened yet)
        std::string fileBaseName = this->PercentEncode(std::string(storableNode->GetName()));
        std::string extension = storageNode->GetDefaultWriteFileExtension();
        std::string storageFileName = fileBaseName + std::string(".") + extension;
        vtkDebugMacro("new file name = " << storageFileName.c_str());
        storageNode->SetFileName(storageFileName.c_str());
    }
    else
    {
        std::vector<std::string> pathComponents;
        vtksys::SystemTools::SplitPath(fileName.c_str(), pathComponents);
        // new file name is encoded to handle : or / characters in the node names
        std::string fileBaseName = this->PercentEncode(pathComponents.back());
        pathComponents.pop_back();
        this->OriginalStorageNodeDirs[storageNode] = vtksys::SystemTools::JoinPath(pathComponents);

        std::string defaultWriteExtension = std::string(".")
            + vtksys::SystemTools::LowerCase(storageNode->GetDefaultWriteFileExtension());
        std::string uniqueFileName = fileBaseName;
        std::string extension = storageNode->GetSupportedFileExtension(fileBaseName.c_str());
        if (defaultWriteExtension != extension)
        {
            // for saving to MRB all nodes will be written in their default format
            uniqueFileName = storageNode->GetFileNameWithoutExtension(fileBaseName.c_str()) + defaultWriteExtension;
        }
        storageNode->SetFileName(uniqueFileName.c_str());
    }

    // also clear out the file list since it's assumed that the default write format is a single file one
    storageNode->ResetFileNameList();
    storageNode->SetDataDirectory(dataDir.c_str());
    vtkDebugMacro("set data directory to "
        << dataDir.c_str() << ", storable node " << storableNode->GetID()
        << " file name is now: " << storageNode->GetFileName());
    // deal with existing files by creating a numeric suffix
    if (vtksys::SystemTools::FileExists(storageNode->GetFileName(), true))
    {
        vtkWarningMacro("file " << storageNode->GetFileName() << " already exists, renaming!");

        std::string uniqueFileName = this->CreateUniqueFileName(fileName);

        vtkDebugMacro("found unique file name " << uniqueFileName.c_str());
        storageNode->SetFileName(uniqueFileName.c_str());
    }

    storageNode->WriteData(storableNode);
}

//----------------------------------------------------------------------------
std::string vtkSlicerPlannerLogic::PercentEncode(std::string s)
{
    std::string validchars = "-_.,@#$%^&()[]{}<>+=";
    std::ostringstream result;

    for (size_t i = 0; i < s.size(); i++)
    {
        if ((s[i] >= 'A' && s[i] <= 'z')
            ||
            (s[i] >= '0'&& s[i] <= '9')
            ||
            (validchars.find(s[i]) != std::string::npos))
        {
            result << s[i];
        }
        else
        {
            result << '%' << std::hex << (unsigned short)s[i] << std::dec;
        }
    }
    return result.str();
}

//----------------------------------------------------------------------------
std::string vtkSlicerPlannerLogic::CreateUniqueFileName(std::string &filename)
{
    std::string uniqueFileName;
    std::string baseName = vtksys::SystemTools::GetFilenameWithoutExtension(filename);
    std::string extension = vtksys::SystemTools::GetFilenameLastExtension(filename);

    bool uniqueName = false;
    int v = 1;
    while (!uniqueName)
    {
        std::stringstream ss;
        ss << v;
        uniqueFileName = baseName + ss.str() + extension;
        if (vtksys::SystemTools::FileExists(uniqueFileName.c_str()) == 0)
        {
            uniqueName = true;
        }
        else
        {
            ++v;
        }
    }
    return uniqueFileName;
}