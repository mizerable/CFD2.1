{-# LANGUAGE ScopedTypeVariables, RankNTypes ,DeriveGeneric , FlexibleInstances #-}
module SolutionDomain where

import qualified Data.Map.Strict as Map
import qualified Data.Set as Set
import Control.Monad
import Data.Maybe
import Data.List
import System.Random.Shuffle
import System.Random
import Data.Array.IO
import qualified Data.Vector as V
import qualified Data.Vector.Mutable as M


import GHC.Generics (Generic)

import Control.DeepSeq


instance Control.DeepSeq.NFData AdjNode
instance Control.DeepSeq.NFData AdjGraph
instance Control.DeepSeq.NFData Position
instance Control.DeepSeq.NFData Grids

data Sign = Positive | Zero | Negative deriving (Enum,Ord,Eq,Show)

data Side =  Now | Prev |East | West | North | South | Top | Bottom | Center deriving (Show,Eq,Ord, Enum)
  
data SchemeType = Directional | Nondirectional       
  
data Direction = Time | X | Y | Z deriving (Enum,Ord,Eq,Show)
data DimensionType = Temporal | Spatial deriving ( Eq)
data Equation a = Equation{
    lhs:: ![a]
    ,rhs:: ![a]
    ,unknownProperty:: !Property}    
data IntegralType = Body | Surface deriving (Eq)
data Property = Speed | Vorticity | Pressure| U | V | W | Density | Temperature | Mew | Dye deriving (Ord,Eq,Enum,Show)
data Position = Position {
    spatialPos:: ![Int]
    ,spatialDimens :: !Int 
    ,timePos:: !Int } deriving (Eq, Ord, Show, Generic)
        
data ValSet a = ValSet{
    calculatedPositions:: ![Position]
    ,vals:: !(Map.Map Position (Map.Map Property a))
    ,areaVal:: !(Map.Map Position (Map.Map Side a))
    ,sideLen:: !(Map.Map Position (Map.Map Direction a))
    ,timeLevelAbsolute :: Int }

data AdjNode  = AdjNode {
    faceArea :: !(V.Vector (Maybe Double ) )--Map.Map Side Double
    ,edgeLen :: !(V.Vector (Maybe Double ) )--Map.Map Direction Double
    ,property :: !(V.Vector (Maybe Double ) )--Map.Map Property Double
    ,neighbors :: !(V.Vector (Maybe Int ) )--Map.Map Side Int 
    ,origPos :: !(Position)
    ,active :: !(Bool)
    ,index :: !(Int)
} deriving (Show, Generic)

data AdjGraph = AdjGraph {
    allNodes:: !Grids
    ,currModTime :: !(Int)
    ,currAbsTime :: !(Int)
} deriving (Show, Generic)

data Grids = Grids {
    grids :: !(V.Vector (V.Vector AdjNode))
} deriving (Show, Generic ) 

sideIndex s = case s of 
    East -> 0
    West ->1
    North ->2
    South ->3
    Top ->4
    Bottom ->5
    Center -> 6
    Now -> 7
    Prev -> 8
    _-> error $ "side index not defined for this " ++ show s
    
dirIndex = getDirectionComponentIndex

propIndex m p = case p of
    U -> 0
    V ->1
    W ->2
    Density ->3
    Temperature ->4
    Mew ->5
    Dye ->6
    _-> error $ "property not stored, maybe it's derived property " ++ show p ++" " ++ m

type AdjPos = (Int, Int) -- ( idx in vector, time position of vector modular time )  

data Expression a = Expression{getTerms:: ![Term a]} 
data Term a = 
    Constant {val:: !a}
    | Unknown { coeff:: !a }
    | SubExpression {expression:: !(Expression a) }  
   -- | Derivative { denom:: !Direction ,function:: !(SchemeType -> Position->Side->Term a), centered:: !Side, 
   --     multiplier:: !(SchemeType->Position->Side-> a) }
    | Derivative_adj { denom_adj:: !Direction ,function_adj:: !(SchemeType -> AdjPos ->Side->Term a), centered_adj:: !Side, 
        multiplier_adj:: !(SchemeType->AdjPos->Side-> a) }  

instance Functor Term where
    fmap f x = case x of
         Constant c -> Constant $ f c
         Unknown u -> Unknown $ f u
         _ -> error "can't fmap other than unknown and constant"

getSign d = case d of   
    0 -> Zero
    _ -> if d > 0 then Positive else Negative

isConvected :: Property -> Bool
isConvected p = case p of
    Density -> False
    Mew -> False
    Pressure -> False
    --Temperature -> False
    _ -> True

isMomentum :: Property -> Bool
isMomentum p = elem p $ enumFrom U \\ enumFrom Density

convectionFromDirection :: Direction -> Property 
convectionFromDirection d = case d of
    X -> U
    Y -> V
    Z -> W
    _ -> error "no convection for that direction"    
 
directionFromConvection d = case d of
    U -> X
    V -> Y
    W -> Z
    _ -> error "no convection for that direction"    

storedSteps:: Int
storedSteps = 3

maxPos :: Direction -> Int
maxPos  d = case d of 
    X -> 100
    Y -> 50
    Z -> 0
    Time -> error "no max time position"
    
gridSize :: Direction -> Double
gridSize d = case d of 
    X -> 1.6
    Y -> 0.8
    Z -> 1
    Time -> error "gridsize for time is the timestep"

cellLength :: Direction -> Double
cellLength d = gridSize d / fromIntegral (maxPos d + 1)

getDirectionComponentIndex :: forall a. Num a => Direction -> a
getDirectionComponentIndex direction = case direction of
    X -> 0
    Y -> 1
    Z -> 2
    _ -> error "not defined or not added dimension"

coordToPos :: [Double] -> Int -> Position
coordToPos coords time = 
    let dirs = enumFrom X
    in Position 
        (map (\ x -> round $ coords !! x / cellLength (dirs !! x))  [0 .. length coords - 1] )
        (length dirs)
        time

posToCoord :: Position -> Direction -> Double
posToCoord (Position s _ _ ) direction = 
    cellLength direction *  fromIntegral (s !! getDirectionComponentIndex direction) 

boundaryPair:: Direction -> (Side,Side)
boundaryPair d = case d of 
     X -> (East,West)
     Y -> (North,South)
     Z -> (Top,Bottom)
     Time ->(Now,Now) -- error "tryng to detect where this is used"

faceToDirections :: Side -> [Direction] 
faceToDirections s = case s of
    East -> [Y,Z]
    West -> [Y,Z]
    North -> [X,Z]
    South -> [X,Z]
    Top -> [X,Y]
    Bottom -> [X,Y]
    Now -> error "no directions outline the now face"
    Prev -> error "no directions outline the prev face"
    Center -> enumFrom X

removeItems  :: (Ord a, Eq a)=> [a] -> [a]-> [a]
removeItems orig remove= 
    let removeSet = Set.fromList remove
    in filter ( `Set.notMember` removeSet) orig      

boundaryPositions :: [Position]
boundaryPositions = obstacles ++ inflowPositions

boundaryPositionsSet :: Set.Set Position
boundaryPositionsSet = Set.fromList $! boundaryPositions 

obstaclesSet :: Set.Set Position
obstaclesSet = Set.fromList $! obstacles  

obstacle :: Position
obstacle = coordToPos [gridSize X / 4 , gridSize Y / 2 , 0 ] 0

obstacle2 :: Position
obstacle2 = Position [quot (maxPos X) 4, quot (maxPos Y) 2 + 5, 0] 3 0

obstacle3 :: Position
obstacle3 = Position [quot (maxPos X) 4, quot (maxPos Y) 2 + 2, 0] 3 0

squareBoundsPts :: [Position]
squareBoundsPts = [
   -- obstacle,
   -- offsetPosition (coordToPos [gridSize X / 4 , gridSize Y / 2 , 0 ] 0) West,
    --coordToPos [gridSize X / 4 - 0.12, gridSize Y / 2  , 0 ] 0
     coordToPos [gridSize X / 4 , gridSize Y / 2 + 0.05 , 0 ] 0
    , coordToPos [gridSize X / 4 + 0.1, gridSize Y / 2 + 0.05 , 0 ] 0
    --, coordToPos [gridSize X / 4 + 0.12, gridSize Y / 2  , 0 ] 0
    , coordToPos [gridSize X / 4 + 0.1, gridSize Y / 2 - 0.05 , 0 ] 0
    , coordToPos [gridSize X / 4 , gridSize Y / 2 - 0.05 , 0 ] 0
    ] 

squareBounds :: [Position] 
squareBounds = connectBounds squareBoundsPts

obstacles :: [Position]
obstacles = 
    let filled =  fillObInterior (Set.fromList squareBounds)  [offsetPosition obstacle East]
        filledGaps = fillEnclosedGaps makeAllPositions $ Set.fromList filled
    in  -- squareBoundsPts
        --squareBounds
        filled ++ filledGaps

timeStep :: Double            
timeStep = 0.00001 --0.0000001

initialMew :: Double
initialMew =  0.1-- 0.000018

initialTemperature :: Double
initialTemperature = 223

sutherlandConstant :: Double
sutherlandConstant = 120

sutherlandLambda :: Double
sutherlandLambda = 
    initialMew *( initialTemperature + sutherlandConstant) / (initialTemperature**1.5)

gasConstantR :: Double
gasConstantR = specificHeatCp - specificHeatCv

specificHeatCp :: Double
specificHeatCp = 1005

specificHeatCv :: Double
specificHeatCv = 716

initialGridPre:: ValSet Double
initialGridPre= 
    let vMap (Position coords _ _) = foldl' (\prev nxt -> Map.insert nxt 
            (case nxt of 
                U-> 1
                V-> 0
                W-> 0
                Density -> 0.3 -- 1.2 
                Mew -> initialMew
                Temperature -> initialTemperature
                Dye-> let midY = div (maxPos Y) 2
                    in case coords of
                        (0 : (y: _ ) ) -> if y == midY then 1 else 0
                        _ -> 0
                _ -> 0 
            ) 
            prev) Map.empty (enumFrom U)
        avMap = foldl' (\prev nxt ->
                    Map.insert nxt (
                        foldl' (\prev2 next2 -> prev2 * fromJust (Map.lookup next2 slMap) ) 1 $ faceToDirections nxt 
                    ) $! prev
                ) Map.empty $!  enumFrom East
        slMap = foldl' (\prev nxt -> 
                    Map.insert nxt (cellLength nxt ) $! prev
                ) Map.empty $! enumFrom X
        v = foldl' (\prev nxt -> Map.insert nxt (vMap nxt) $! prev) Map.empty  $!  makeAllPositions
        av = foldl' (\prev nxt -> Map.insert nxt avMap $! prev) Map.empty $! makeAllPositions
        sl = foldl' (\prev nxt -> Map.insert nxt slMap $! prev) Map.empty $! makeAllPositions
    in ValSet makeAllPositions v av sl 0 
  
initialGrid :: ValSet Double    
initialGrid = 
    let calcPos = removeItems (calculatedPositions initialGridPre) $! boundaryPositions
        v= vals initialGridPre
        av = areaVal initialGridPre
        sl = sideLen initialGridPre
    in foldl'
        (\prev nextPos ->
            foldl'
                (\prev2 nextProp -> setVal prev2 nextPos nextProp 0.0) 
                prev
                $ enumFrom U \\ enumFrom Density     
        )
        (ValSet calcPos v av sl 0)
         $! obstacles

pseudoRandList size = 
    let (list, _) = foldl' (\(l,s) lim ->
                let nr =randLessThan lim s  
                in ( l++[div lim 2],nr) 
            ) ([],1) $ reverse [1..size]
    in list
    
randLessThan lim n = 2 * n `mod` lim  

initialGrid_adj :: AdjGraph
initialGrid_adj = 
    let cp = Set.fromList $ calculatedPositions initialGrid
        toVector im d =
            let inputs = Map.toList d
            in V.foldl' 
                (\prev (k,a)->prev V.// [( im k , Just a)] ) 
                (V.replicate (length inputs) Nothing) 
                $ V.fromList inputs
        (list_nodes_unconnected,locations_list,_) = foldl'
            (\(list,m,i) next_pos ->
                (list ++ 
                    [AdjNode
                        (toVector sideIndex $ fromJust $ Map.lookup next_pos $ areaVal initialGrid) --face area
                        (toVector dirIndex $ fromJust $ Map.lookup next_pos $ sideLen initialGrid) -- side length
                        (toVector (propIndex "from initializing") $ fromJust $ Map.lookup next_pos $ vals initialGrid) -- values
                        (V.replicate (length $ enumFrom East) Nothing)--neighbors
                        next_pos --origpos
                        (Set.member next_pos cp)--active
                        i --index
                    ]    
                , m++[(next_pos,i)] 
                , i+1)
            )
            ([],[],0) makeAllPositions
        locations = foldl' (\prev (key,val)->Map.insert key val prev) Map.empty 
            $ shuffle locations_list $ pseudoRandList $ length locations_list -1
        list_nodes = 
            map 
                (\(AdjNode a b c neighbs op e i) ->
                    AdjNode 
                        a b c
                        (
                            foldl'
                                (\prev x -> 
                                    let n = offsetPosition op x
                                        n_idx = Map.lookup n locations 
                                    in case n_idx of 
                                            Just l -> prev V.// [(sideIndex x , Just l )] 
                                            Nothing ->error $ "can't find previously added location when translating valset to adjgraph: " 
                                                ++ show op ++ " | neighbor: " ++ show n ++ " | locations: " ++ show locations
                                    --in Map.insert x 888 prev
                                ) neighbs $ enumFrom East
                        ) 
                        op e i
                ) list_nodes_unconnected   
        aln = V.replicate storedSteps (V.fromList list_nodes) 
    in AdjGraph (Grids $ force aln) 0 0
    
setVal:: ValSet Double -> Position -> Property -> Double -> ValSet Double
setVal (ValSet p v av sl tl) pos prp newVal = 
    let subDict = case Map.lookup pos v  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p (Map.insert pos (Map.insert prp newVal subDict) v) av sl tl

setFaceArea :: ValSet Double -> Position -> Side -> Double -> ValSet Double
setFaceArea (ValSet p v av sl tl) pos side newVal = 
    let subDict = case Map.lookup pos av  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p v (Map.insert pos (Map.insert side newVal subDict) av) sl tl             
    
setCellLength :: ValSet Double -> Position -> Direction -> Double -> ValSet Double
setCellLength (ValSet p v av sl tl) pos dir newVal = 
    let subDict = case Map.lookup pos sl  of
            Nothing -> Map.empty
            Just sb -> sb
    in ValSet p v av (Map.insert pos (Map.insert dir newVal subDict) sl) tl
        
cartProd:: [[a]] -> [[a]] -> [[a]]
cartProd xs ys = [ x ++ y | x <- xs, y <- ys]

inflowPositions::[Position]
inflowPositions =  makePositions $!  0 :  (map maxPos  $! enumFrom Y) 

makeAllPositions::[Position]
makeAllPositions = makePositions $! map maxPos $! enumFrom X 
         
makeAllPositions_adj :: [AdjPos]         
makeAllPositions_adj = map (\x-> (x,0)) [0.. length makeAllPositions - 1]         
         
makePositions :: [Int] -> [Position]
makePositions maxes =
    let ranges = map (\x -> map (:[]) [0..x] ) maxes 
        posCoords = foldl' cartProd [[]] ranges
    in map (\coords -> Position coords (length coords) 0) posCoords

setElem newElem idx list = 
    map (\x -> if x == idx then newElem else list!!x) [0..length list -1]

modifyPositionComponent::Position -> Direction -> Int -> Position
modifyPositionComponent (Position p d t) direction amt= case direction of 
    Time -> Position p d amt
    _ -> Position (setElem amt (getDirectionComponentIndex direction) p) d t

isUpperSide:: Side -> Bool
isUpperSide side = case side of
    East -> True
    North -> True
    Top -> True
    Now -> True 
    _->False

directionFromCenter:: Side-> Direction
directionFromCenter side = case side of
    East -> X
    West -> X
    North -> Y
    South -> Y
    Top -> Z
    Bottom -> Z
    Now -> Time
    Prev -> Time
    Center -> error "there is no direction from center" 

pushBackTime :: Int -> Int
pushBackTime t = mod (t - 1) storedSteps

pushUpTime :: Int -> Int
pushUpTime t = mod (t+1) storedSteps    
    
advanceTime :: Position -> Position
advanceTime (Position p d t ) = Position p d $! pushUpTime t      

advanceTime_adj :: AdjPos -> AdjPos
advanceTime_adj (idx,t) = (idx, pushUpTime t)
        
offsetPosition:: Position->Side ->Position
offsetPosition (Position p d t) side = case side of
    Center -> Position p d t 
    Now -> Position p d t 
    Prev -> Position p d $! pushBackTime t   
    _ -> 
        let position = Position p d t
            maxOrMin = if isUpperSide side then min else max
            offsetAmount = if isUpperSide side then 1 else (-1)
            direction = directionFromCenter side
            boundary = if isUpperSide side 
                then maxPos direction
                else 0
        in modifyPositionComponent position direction 
            $ maxOrMin boundary $ getPositionComponent position direction + offsetAmount   

getNode :: AdjGraph -> AdjPos -> AdjNode
getNode (AdjGraph g _ _) (i,t) = (grids g) V.! t V.! i

-- used by central differencing,. ONLY FOR NON DERIVED 
getNodeVal :: AdjGraph -> AdjPos -> Property -> Maybe Double
getNodeVal g x p =
    let v = property $ getNode g x
    in v V.! (propIndex "from getNodeVal(only in central diff)" p)  

offsetPosition_adj:: AdjPos->Side -> AdjGraph ->AdjPos
offsetPosition_adj (idx,t) side env = case side of
    Center -> (idx,t) 
    Now -> (idx,t) 
    Prev -> error "can't push back to another time level" --(idx, pushBackTime t)   
    _ -> 
        let neighs = (neighbors $ getNode env (idx,t)) V.! (sideIndex side)
        in case neighs of
            Just n -> (n,t)
            Nothing -> (idx,t)

getPositionComponent:: Position -> Direction -> Int
getPositionComponent (Position p _ t) direction = case direction of 
    Time -> t
    _ -> p !! getDirectionComponentIndex direction

average::(Num a, Fractional a) => [a]->a
average terms =
    let len = fromIntegral $! length terms 
        f p n= p + n / len
    in foldl' f 0 terms

isBoundaryPosition :: Position -> Bool
isBoundaryPosition (Position p d _) = Set.member (Position p d 0) boundaryPositionsSet

isObstaclePosition :: Position -> Bool
isObstaclePosition (Position p d _) = Set.member (Position p d 0) obstaclesSet 

positionIfWall :: Position -> Position
positionIfWall (Position p d t) = if isBoundaryPosition (Position p d t) 
    then Position p d 0
    else Position p d t
    
envIfWall :: Position -> ValSet Double -> ValSet Double
envIfWall (Position p d _) env = if isBoundaryPosition (Position p d 0) 
    then initialGrid
    else env       

pecletNumber :: Position -> Side -> ValSet Double -> Double
pecletNumber position side env = 
    let direc = directionFromCenter side
        momentum = propCentralDiff (convectionFromDirection direc) position side env  
        density = propCentralDiff Density position side env
        viscosity = propCentralDiff Mew position side env
        l = sideLength direc position env 
    in (density * momentum * l) / viscosity 

approximateDerivative_adj :: Term Double -> AdjPos -> AdjGraph -> Term Double
approximateDerivative_adj deriv position vs= case deriv of 
    (Derivative_adj direction func side m) ->
        if directionFromCenter side == direction 
        then
            let neighbor = offsetPosition_adj position side vs
                sl =  sideLength_adj direction position vs 
                sln = sideLength_adj direction neighbor vs
                interval = average [sl, sln ]
                thisVal = func Nondirectional position Center 
                neighborVal = func Nondirectional neighbor Center
                neighborIsUpper = isUpperSide side  
                f first = if first==neighborIsUpper then neighborVal else thisVal
                mult = m Nondirectional position Center
            in case (f True, f False) of
                (Constant c1 , Constant c2) -> Constant $ (c1-c2)*mult/interval 
                _ -> error "can't approx >1 order derivs. deriv must produce constants" 
        else
            let (s1, s2) = boundaryPair direction
                n1 = offsetPosition_adj position s1 vs
                n2 = offsetPosition_adj position s2 vs
                sl =  sideLength_adj direction position vs
                sln1 = sideLength_adj direction n1 vs
                sln2 = sideLength_adj direction n2 vs
                interval = 2 * average [sl,sln1,sln2]
                n1Val = func Nondirectional n1 side
                n2Val = func Nondirectional n2 side
                mult = m Nondirectional position side
            in case (n1Val, n2Val) of
                (Constant c1 , Constant c2) -> 
                    Constant $ (c1-c2)*mult/interval  
                _ -> error "can't approx >1 order derivs. deriv must produce constants"
    _ -> error "can't approx something that isn't a deriv"

magnitude :: forall s. Floating s => [s] -> s
magnitude v = sqrt $ foldl' (\prev nVal ->prev + (nVal*nVal)) 0.0 v

prop_adj :: SchemeType -> Property -> AdjPos -> Side -> AdjGraph -> Double
prop_adj schemeType =
    let f scheme prp pos side env  = 
            if side == Center || side == Now || side == Prev 
                then propCentralDiff_adj prp pos side env 
                else scheme prp pos side env 
    in 
        let scheme = case schemeType of
                Directional -> f propLimitedSlope_adj
                Nondirectional -> f propLimitedSlope_adj
            momentums = enumFrom U \\ enumFrom Density
        in \prp pos side env-> case prp of
                Speed-> magnitude $ map (\nxt -> scheme nxt pos side env) momentums   
                Vorticity -> 
                    let pairs = [(x,y) | x <- momentums , y <- tail $ dropWhile (/= x) momentums ]
                        vortComponents = map (\(a,b) ->
                                let makeDeriv m n 
                                        = Derivative_adj (directionFromConvection m) 
                                            (\_ p s -> Constant $ scheme n p s env) East -- east is arbitrary here  
                                            (\_ _ _ -> 1) 
                                    deriv1 = makeDeriv a b
                                    deriv2 = makeDeriv b a
                                in val (approximateDerivative_adj deriv1 pos env) -
                                    val (approximateDerivative_adj deriv2 pos env)    
                            ) pairs
                    in magnitude vortComponents 
                Mew -> 
                    let t = scheme Temperature pos side env 
                    in sutherlandLambda * (t**1.5)/ (t + sutherlandConstant)
                Pressure ->
                    gasConstantR * scheme Temperature pos side env *
                        scheme Density pos side env
                _-> scheme prp pos side env  

propDirectional prp position side env =
    let neighbor = offsetPosition position side
        decide =
            let peclet = pecletNumber position side env  
            in (if (peclet == 0) || (abs peclet < 0.99) 
                then propCentralDiff 
                else(if abs peclet > 2.6 
                        then propUpwindDiff peclet 
                        else propQUICK peclet)) 
    in case (isObstaclePosition neighbor || isObstaclePosition position
                , isMomentum prp 
                ,elem side ( enumFrom East \\ enumFrom Center )
                    && isConvected prp ) of
        --(True,True,_)-> 0.0
        (_,_,True) -> decide prp position side env 
        _ -> propCentralDiff prp position side env 
        
oppositeSide :: Side -> Side 
oppositeSide s = 
    let d = directionFromCenter s
        (s1,s2) = boundaryPair d
    in if isUpperSide s 
        then s2 else s1 

propQUICK :: Double -> Property -> Position -> Side ->ValSet Double -> Double
propQUICK pec prp position side env = 
    let upper = if isUpperSide side then side else Center 
        lower = if isUpperSide side then Center else side
        doubleOffset pos = offsetPosition pos >>= offsetPosition
        farPoint upBracket downBracket = if upBracket == Center
            then offsetPosition position $ oppositeSide downBracket
            else doubleOffset position upBracket 
    in if pec >= 0
        then ((6/8)* propCentralDiff prp (offsetPosition position lower) Center env ) 
            + ((3/8)* propCentralDiff prp (offsetPosition position upper) Center env ) 
            - ((1/8)* propCentralDiff prp (farPoint lower upper) Center env ) 
        else ((6/8)* propCentralDiff prp (offsetPosition position upper) Center env ) 
            + ((3/8)* propCentralDiff prp (offsetPosition position lower) Center env ) 
            - ((1/8)* propCentralDiff prp (farPoint upper lower) Center env )

propUpwindDiff :: Double -> Property -> Position -> Side -> ValSet Double -> Double              
propUpwindDiff pec prp position side env = 
    let upper = if isUpperSide side then side else Center 
        lower = if isUpperSide side then Center else side
    in if pec >= 0
        then propCentralDiff prp (offsetPosition position lower) Center env
        else propCentralDiff prp (offsetPosition position upper) Center env

-- http://www.ammar-hakim.org/_static/files/1010-muscl-hancock.pdf
-- http://bulletin.pan.pl/(60-1)45.pdf
propLimitedSlope_adj prp position side env = 
    let valCenter = propCentralDiff_adj prp position Center env
        d = directionFromCenter side
        interval = sideLength_adj d position env
        (upper, lower) = boundaryPair d
        upperNVal:(lowerNVal:_)  
            = map (\s -> propCentralDiff_adj prp (offsetPosition_adj position s env) Center env) 
                [upper,lower] 
        ave = superbee (upperNVal - valCenter) (valCenter - lowerNVal)
        --ave = epsilon (upperNVal - valCenter) (valCenter - lowerNVal) (interval *interval *interval)
        --ave = minmodLimit (upperNVal - valCenter) (valCenter - lowerNVal)
        -- ave = ospre (upperNVal - valCenter) (valCenter - lowerNVal)
        -- ave = vanLeer (upperNVal - valCenter) (valCenter - lowerNVal)
        --ave = ((upperNVal - valCenter)+(valCenter - lowerNVal))/2
    in (if isUpperSide side then (+) else (-)) valCenter  
            $ 0.5 * ave   

propCentralDiff :: Property -> Position -> Side -> ValSet Double -> Double
propCentralDiff prp position side env = 
    let neighbor = offsetPosition position side
        noValError = error ("no value "
                            ++ 
                                foldl' (\prev nxt -> prev ++ " " ++ show (spatialPos position !! nxt)) "" 
                                    [0..spatialDimens position-1] 
                            ++ " "
                            ++ show (timePos position)++ " "
                            ++ show prp ++ " "
                            ++ show side)
        getVal:: Position -> Map.Map Position (Map.Map Property Double) -> Double
        getVal p set = fromMaybe 
            --noValError
            ( case  Map.lookup (modifyPositionComponent p Time 0) (vals initialGrid )>>= Map.lookup prp of
                        Nothing -> noValError
                        Just r -> r
            )   
            (Map.lookup p set >>= Map.lookup prp)
        res p = getVal p (vals env )
    in case (isObstaclePosition neighbor || isObstaclePosition position
                , isMomentum prp, position == neighbor) of
        --(True,True,_) -> 0.0
        (_,_,True) -> res position
        _ -> average [ res position, res neighbor]

propCentralDiff_adj :: Property -> AdjPos -> Side -> AdjGraph-> Double
propCentralDiff_adj prp (idx,t) side env = 
    let neigh = offsetPosition_adj (idx,t) side env
        neighborPos = origPos $ getNode env neigh
        originalPos = origPos $ getNode env (idx,t)
        noValError = error ("no value "
                            ++ 
                                foldl' (\prev nxt -> prev ++ " " ++ show (spatialPos originalPos !! nxt)) "" 
                                    [0..spatialDimens originalPos-1] 
                            ++ " "
                            ++ show t ++ " "
                            ++ show prp ++ " "
                            ++ show side)
        getVal (i_1,t_1) = case getNodeVal env (i_1,t_1) prp of
            Just s -> s
            Nothing -> case prp of
                Density -> noValError
                Temperature -> noValError
                Pressure -> noValError 
                _ -> case getNodeVal env (i_1, 0) prp of
                    Just s1 -> s1
                    Nothing -> noValError
    in case (isObstaclePosition neighborPos || isObstaclePosition originalPos
                , isMomentum prp, originalPos == neighborPos) of
        --(True,True,_) -> 0.0
        (_,_,True) -> getVal (idx,t)
        _ -> average [ getVal (idx,t), getVal neigh]

-- GEOMETRY STUFF 
 
fillObInterior :: Set.Set Position -> [Position] -> [Position]
fillObInterior obBounds points = case points of
    [] -> Set.toList obBounds
    x:xs -> 
        let unvN = filter (`Set.notMember` obBounds) $ getNeighbors x
        in fillObInterior (Set.union obBounds $ Set.fromList unvN) $ xs ++ unvN
        
distanceSquared :: Position -> Position -> Double        
distanceSquared (Position sp1 _ _) (Position sp2 _ _) = 
    let components = zip sp1 sp2         
    in fromIntegral $ foldl' (\prev (a,b) -> prev + (a-b)*(a-b) ) 0 components
    
distance:: Position -> Position -> Double    
distance p1 p2 = sqrt $ distanceSquared p1 p2 
        
isBetween :: Position-> Position -> Position -> ValSet Double -> Bool
isBetween p1 p2 testPt vs= 
    let aSq = distanceSquared p1 testPt
        bSq = distanceSquared p2 testPt
        cSq = distanceSquared p1 p2
        x = sqrt $ (aSq + bSq - cSq ) / 2
        cornerDist = sqrt $ 
            foldl' (\prev nxt ->
                    let sl = sideLength nxt testPt vs
                    in prev + (sl*sl) 
                ) 0 (enumFrom X)
    in cornerDist > x
        && testPt /= p1 && testPt /= p2        
        
connectTwo p1 p2 allPos vs = filter (isBetween p1 p2 vs) allPos

getNeighborsCorners :: Position -> [Position]
getNeighborsCorners position =
    let cross = filter (uncurry (/=) ) [(x,y) | x <- enumFrom East, y<-enumFrom East] 
        n = nub $ filter (/= position)
            $ map (\(a,b) -> offsetPosition (offsetPosition position a) b ) cross
    in shuffle' n (length n) $ mkStdGen $ sum $ spatialPos position

getNeighbors :: Position -> [Position]
getNeighbors position =
    nub $ filter (/= position)
        $ map (offsetPosition position) $ enumFrom East

tracePath :: forall t a.
               Ord a =>
               Map.Map a (t, a) -> a -> [a] -> [a]
tracePath g end path = 
    let prev = Map.lookup end g
    in case prev of
        Nothing -> path
        Just (_, p) -> if p == end 
            then path 
            else tracePath g p $ p:path  

hugeNumber :: Double
hugeNumber = 99999999999.9

rounder :: forall a s b.(RealFrac s, Integral b, Fractional a) => s -> b -> a
rounder f n = fromInteger (round $ f * (10 ^ n)) / (10.0 ^^ n)

shortestPath :: Map.Map Position (Double, Position)
                  -> Position -> Map.Map Position (Double, Position) -> [Position]
shortestPath unvisited end visited = if unvisited==Map.empty 
    then error "problem in shortest path.. searched everything and wasn't found" 
    else 
        let comparer a= if mod (Map.size unvisited * Map.size visited + a ) 2 < 1 
                then (<) else (<=)
            (pos, (dist,predecessor )) = 
                foldl' (\(p1,(p2,p3) ) (n1,(n2,n3)) -> if comparer (fromIntegral $ round $ n2*p2*57) n2  p2 
                            then (n1,(n2,n3)) else (p1,(p2,p3))  ) 
                        (end,(hugeNumber,end)) 
                        $ shuffle' (Map.toList unvisited) (Map.size unvisited) $ mkStdGen 137
            neighbs = filter (`Map.notMember` visited) $ getNeighborsCorners pos 
            updatedUnvisited = foldl' (\prev nxt -> 
                    let newd = distance nxt pos + dist
                        oldd = fst $ fromMaybe (hugeNumber,end) (Map.lookup nxt prev) 
                    in if comparer (fromIntegral $ round $ newd*oldd*57) newd oldd 
                        then Map.insert nxt (newd,pos) prev
                        else prev
                ) (Map.delete pos unvisited) neighbs
            updatedVisited = Map.insert pos (dist,predecessor ) visited
        in (if pos == end 
            then tracePath updatedVisited end [] 
            else shortestPath updatedUnvisited end updatedVisited)
           
connectBounds :: [Position] -> [Position]
connectBounds points = fst $ foldl' 
    (\(prevPts, lastPt) nxt -> 
        ( prevPts 
            ++ shortestPath 
                (Map.insert lastPt (0.0,lastPt) Map.empty) 
                nxt 
                Map.empty 
        , nxt) ) 
    ([],last points) points

fillEnclosedGaps :: [Position] -> Set.Set Position -> [Position]
fillEnclosedGaps allPos wallPos = 
    let added = filter
            (\x ->
                let opposingWalls = foldl'
                        (\prev nxt->
                            let (s1 ,s2) = boundaryPair nxt
                            in prev ||
                                ( Set.member (offsetPosition x s1) wallPos 
                                    && Set.member (offsetPosition x s2) wallPos )
                        )
                        False
                        $ enumFrom X
                in Set.notMember x wallPos && opposingWalls
            )
            allPos
    in case added of 
        [] -> []
        _ -> added ++ fillEnclosedGaps allPos (Set.union wallPos $ Set.fromList added) 

sideArea :: forall a. Num a => Side -> Position -> ValSet a -> a
sideArea s (Position p d _) vs = case s of 
    Now -> 1
    Prev -> 1
    _ -> fromMaybe  (error $ "error getting side area for " ++ show p ++ " " ++ show s)
            (Map.lookup (Position p d 0) (areaVal $! vs) >>= Map.lookup s)

sideArea_adj :: Side -> AdjPos -> AdjGraph -> Double
sideArea_adj side x env = case side of 
    Now -> 1
    Prev -> 1
    _ ->
        let node = getNode env x
            v = faceArea node
            fa = v V.! (sideIndex side)
        in case fa of
            Just a -> a
            Nothing -> error $ "error getting side area for " ++ show (origPos node) ++ " " ++ show side

sideLength :: Direction -> Position -> ValSet Double -> Double
sideLength direction (Position p d _) vs = case direction of 
    Time -> timeStep
    _ -> fromMaybe  (error $ "error getting side length for " ++ show p ++ " " ++ show direction)
            (Map.lookup (Position p d 0) (sideLen $! vs) >>= Map.lookup direction)    
            
sideLength_adj :: Direction -> AdjPos -> AdjGraph -> Double
sideLength_adj direction x env = case direction of 
    Time -> timeStep
    _ -> 
        let node = getNode env x
            v = edgeLen node
            sl = v V.! (dirIndex direction) 
        in case sl of
            Just l -> l
            Nothing ->  error $ "error getting side length for " ++ show (origPos node) ++ " " ++ show direction

chooseSlopeHelper :: forall a a1 s.
                       (Ord a1, Ord a, Num a1, Num a) =>
                       (a -> a1 -> s) -> (a -> a1 -> s) -> a -> a1 -> Maybe s
chooseSlopeHelper f1 f2 x y =
    let sign = getSign x 
    in if getSign x == getSign y
        then case sign of
                Positive -> Just $ f1 x y
                Negative -> Just $ f2 x y
                Zero -> Nothing 
        else Nothing 

chooseSlope :: forall b.
                 (Ord b, Fractional b) =>
                 (b -> b -> b) -> (b -> b -> b) -> [b] -> b
chooseSlope f1 f2 n = 
    let res = foldM
            (chooseSlopeHelper f1 f2) 
            (head n) n
    in fromMaybe 0.0 res 

minmod :: [Double] -> Double
minmod = chooseSlope min max

maxmod :: [Double] -> Double
maxmod = chooseSlope max min

superbee :: Double -> Double -> Double
superbee a b = minmod [maxmod [a,b], minmod [2*a,2*b] ]

minmodLimit :: Double -> Double -> Double
minmodLimit a b = minmod [ (a + b) /2 , 2*a, 2*b ]

epsilon :: forall a. Fractional a => a -> a -> a -> a
epsilon a b eSq = ( ((b*b + eSq )*a)  + ((a*a+eSq)*b) ) / ((a*a) + (b*b) + (2*eSq)) 

vanLeer :: Double -> Double -> Double
vanLeer = altFormLimiter (\r -> let abs_r = abs r in (r+abs_r)/(1+abs_r) )  

ospre :: Double -> Double -> Double
ospre = altFormLimiter (\r -> 1.5 * (r*r + r) / (r*r + r + 1) )

altFormLimiter :: forall a.
                    (Ord a, Fractional a) =>
                    (a -> a) -> a -> a -> a
altFormLimiter rFunc a b = 
    let expr steep shallow = 
            let r = steep/shallow
            in shallow * rFunc r
    in case (getSign a, getSign b) of
        (Positive,Positive) -> expr (max a b) (min a b)
        (Negative,Negative) -> expr  (min a b) (max a b)
        _-> 0.0 


       