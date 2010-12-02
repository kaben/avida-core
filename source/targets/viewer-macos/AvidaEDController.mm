//
//  AvidaEDController.mm
//  Avida
//
//  Created by David on 11/30/10.
//  Copyright 2010 Michigan State University. All rights reserved.
//
//
//  This file is part of Avida.
//
//  Avida is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
//  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
//  Avida is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License along with Avida.
//  If not, see <http://www.gnu.org/licenses/>.
//

#import "AvidaEDController.h"

#import "AvidaAppDelegate.h"
#import "AvidaRun.h"
#import "MapGridView.h"


@implementation AvidaEDController

- (id)initWithAppDelegate: (AvidaAppDelegate*)delegate {
  self = [super initWithWindowNibName:@"Avida-ED-MainWindow"];
  
  if (self != nil) {
    app = delegate;

    currentRun = nil;
    listener = NULL;
    
    [self showWindow:self];
  }
  
  return self;
}


- (void) dealloc {
  delete listener;
  listener = NULL;
  [super dealloc];
}


- (void) finalize {
  delete listener;
  listener = NULL;
  [super finalize];
}


- (void) windowDidLoad {
  [btnRunState setState:NSOffState];
}


- (IBAction) toggleRunState:(id)sender {
  if ([sender state] == NSOnState && currentRun == nil) {
    NSFileManager* fileManager = [NSFileManager defaultManager];
    NSArray* urls = [fileManager URLsForDirectory:NSDocumentDirectory inDomains:NSUserDomainMask];
    
    if ([urls count] > 0) {
      NSURL* userDocumentsURL = [urls objectAtIndex:0];
      NSURL* runURL = [NSURL URLWithString:@"../test" relativeToURL:userDocumentsURL];
      currentRun = [[AvidaRun alloc] initWithDirectory:runURL];
    }
    if (currentRun == nil) {
      NSAlert* alert = [[NSAlert alloc] init];
      [alert addButtonWithTitle:@"OK"];
      NSString* msgText = [NSString stringWithFormat:@"Unable to load run configuration in \"test\""];
      [alert setMessageText:msgText];
      [alert setInformativeText:@"Check the run log for details on how to correct your configuration files."];
      [alert setAlertStyle:NSWarningAlertStyle];
      [alert beginSheetModalForWindow:[sender window] modalDelegate:nil didEndSelector:nil contextInfo:nil];
      [sender setState:NSOffState];
    } else {
      if (!listener) listener = new cMainThreadListener(self);
      [currentRun attachListener:self];
      
      [txtUpdate setStringValue:@"Time (updates): 0"];
    }
  } else {
    if ([sender state] == NSOffState) {
      [currentRun resume];
    } else {
      [currentRun pause];
    }
  }
}


- (BOOL) splitView:(NSSplitView*)splitView canCollapseSubview:(NSView*)subview {
  if (splitView == mainSplitView && subview == [[splitView subviews] objectAtIndex:0]) {
    return YES;
  }
  
  return NO;
}


-(CGFloat) splitView:(NSSplitView*)splitView constrainMaxCoordinate:(CGFloat)proposedMax ofSubviewAt:(NSInteger)dividerIndex {
  return proposedMax;
}


-(CGFloat) splitView:(NSSplitView*)splitView constrainMinCoordinate:(CGFloat)proposedMin ofSubviewAt:(NSInteger)dividerIndex {
  if (splitView == mainSplitView && dividerIndex == 0) {
    NSView* subview = [[splitView subviews] objectAtIndex:dividerIndex];
    NSRect subviewFrame = subview.frame;
    CGFloat frameOrigin = subviewFrame.origin.x;
    
    return frameOrigin + 150;
  }
  
  return proposedMin;
}

- (void) windowWillClose: (NSNotification*)notification {
  if (currentRun != nil) {
    [currentRun end];
    currentRun = nil;
  }
  [app removeWindow:self];
}


@synthesize listener;


- (void) handleMap: (CoreViewMap*)pkg {
  [mapView updateState: [pkg map]];
}


- (void) handleUpdate: (CoreViewUpdate*)pkg {
  NSString* str = [NSString stringWithFormat:@"Update: %d", [pkg update]];
  [txtUpdate setStringValue:str]; 
}


@end